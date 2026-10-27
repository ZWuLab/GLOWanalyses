# ==============================================================================
# File: cloud/terra/task-monitor.sh
# Purpose: a resource monitor for the tasks of glow_scan.wdl, sourced by each task's command
#   from the image (/opt/GLOWanalyses/cloud/terra/task-monitor.sh). glow_monitor_start FILE
#   starts a background loop that samples the memory in use every 15 seconds and writes, every
#   10 minutes and once more at exit, one line with the memory in use, the peak so far, the load
#   average and the vCPUs, to FILE and to the task's standard output. The peak is the largest
#   sample, so a task shorter than the interval reports the memory before its work started. The exit line is also
#   written when the command fails, so an out-of-memory stop leaves its last readings in the
#   task log. It is what sizes mem_gb_scan per task kind: the peak of a representative
#   chromosome, with a margin, is the number to request.
#
#   Memory in use is the memory of this task's own processes: the proportional set size (PSS) of
#   the task shell and all its descendants, summed from /proc/<pid>/smaps_rollup, so the pages
#   that forked workers share are counted once; a process without a readable smaps_rollup
#   contributes its resident set instead. It needs no cgroup, which matters because a container
#   runtime without a cgroup namespace (Apptainer on a cluster) shows the whole machine's cgroup
#   at /sys/fs/cgroup. The exit line also prints the machine's total memory, which on Terra is
#   the task's own virtual machine.
#
#   Usage, in a bash command that runs under set -euo pipefail:
#     . /opt/GLOWanalyses/cloud/terra/task-monitor.sh
#     glow_monitor_start "$PWD/resources.txt"     # [interval_s = 15] [lines_every = 40]
#   It installs an EXIT trap that stops the loop and writes the exit line; a task that sets its
#   own EXIT trap afterwards replaces it. The loop is a background job of the task shell, so a
#   bare `wait` in the task would wait for it too: wait for a specific pid instead.
#
#   Authors: Zheyang Wu (Worcester Polytechnic Institute), for ZWuLab. Derived from the
#   monitor loop of the Terra workflows written by Alberto Brusati (Landers lab, UMass Chan
#   Medical School), 2026, who used its readings to size the machines per task kind.

# The memory of this task's process tree, in MB: the PSS of the shell GLOW_MON_ROOT and every
# descendant (RSS when a process's smaps_rollup is not readable). The tree is walked from
# /proc/<pid>/stat (the parent is the field after the command's closing parenthesis), in the
# numbering of the proc mount, which is also GLOW_MON_ROOT's.
glow_mem_mb() {
  local root="${GLOW_MON_ROOT:-$$}" d pid stat rest kb total=0 p c
  local -A parent=() seen=()
  for d in /proc/[0-9]*; do
    pid=${d#/proc/}
    read -r stat 2>/dev/null < "$d/stat" || continue   # the process may have exited since the glob
    rest=${stat##*) }                 # "state ppid ...": the command name may hold spaces and parentheses
    set -- $rest; parent[$pid]=${2:-0}
  done
  local -a queue=("$root"); seen[$root]=1
  while [ ${#queue[@]} -gt 0 ]; do
    p=${queue[0]}; queue=("${queue[@]:1}")
    for c in "${!parent[@]}"; do
      if [ "${parent[$c]}" = "$p" ] && [ -z "${seen[$c]:-}" ]; then seen[$c]=1; queue+=("$c"); fi
    done
  done
  for p in "${!seen[@]}"; do
    kb=$(awk '/^Pss:/ { s += $2 } END { print s + 0 }' "/proc/$p/smaps_rollup" 2>/dev/null || echo 0)
    if [ "${kb:-0}" -gt 0 ]; then total=$(( total + kb ))
    elif read -r _ kb _ 2>/dev/null < "/proc/$p/statm"; then total=$(( total + kb * 4 )); fi   # pages of 4 KB
  done
  echo $(( total / 1024 ))
}

# The machine's total memory, in GB (on Terra, the task's own virtual machine).
glow_mem_total_gb() { awk '/^MemTotal:/ { printf "%.1f", $2 / 1048576 }' /proc/meminfo; }

# One reading as a line: the time, the memory in use and the peak so far (GB), the load average, the vCPUs.
glow_monitor_line() {   # glow_monitor_line <label> <used_mb> <peak_mb>
  printf '[monitor%s %s] memory in use %s GB, peak %s GB, load average %s, vCPUs %s\n' \
    "${1:+ $1}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "$(awk -v m="$2" 'BEGIN { printf "%.1f", m / 1024 }')" "$(awk -v m="$3" 'BEGIN { printf "%.1f", m / 1024 }')" \
    "$(cut -d' ' -f1-3 /proc/loadavg)" "$(nproc)"
}

# Start the loop. The peak crosses from the loop to the exit line through a file, written whole
# (write, then rename), so the exit line never reads a half-written value.
glow_monitor_start() {   # glow_monitor_start <file> [interval_s=15] [lines_every=40]
  GLOW_MON_FILE="$1"; GLOW_MON_PEAK="$1.peak"
  local interval="${2:-15}" every="${3:-40}" self
  # The root of the tree measured is this shell, numbered as /proc numbers it: a builtin read of
  # /proc/self/stat gives that number without a fork, where $$ would give the shell's number in its
  # own PID namespace, which differs from /proc's when a container runtime starts the task with a
  # new PID namespace over the host's proc mount (miniwdl's Apptainer tasks do).
  if read -r self _ 2>/dev/null < /proc/self/stat; then GLOW_MON_ROOT=$self; else GLOW_MON_ROOT=$$; fi
  : > "$GLOW_MON_FILE"; echo 0 > "$GLOW_MON_PEAK"
  (
    set +e
    local_sleep=""
    trap 'kill "$local_sleep" 2>/dev/null; exit 0' TERM   # sleep is a child: on TERM take it along
    peak=0; n=0
    while :; do
      used=$(glow_mem_mb)
      if [ "$used" -gt "$peak" ]; then
        peak=$used; echo "$peak" > "$GLOW_MON_PEAK.tmp" && mv -f "$GLOW_MON_PEAK.tmp" "$GLOW_MON_PEAK"
      fi
      if [ $(( n % every )) -eq 0 ]; then glow_monitor_line "" "$used" "$peak" | tee -a "$GLOW_MON_FILE"; fi
      n=$(( n + 1 ))
      sleep "$interval" & local_sleep=$!
      wait "$local_sleep"
    done
  ) &
  GLOW_MON_PID=$!
  trap glow_monitor_stop EXIT
}

# Stop the loop and write the exit line: the peak, the machine's total, the vCPUs and the command's exit status.
glow_monitor_stop() {
  local rc=$? peak
  kill "$GLOW_MON_PID" 2>/dev/null || true
  wait "$GLOW_MON_PID" 2>/dev/null || true
  peak=$(cat "$GLOW_MON_PEAK" 2>/dev/null || echo 0)
  printf '[monitor exit %s] peak memory in use %s GB (this task'"'"'s processes), machine %s GB, vCPUs %s, exit status %d\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(awk -v m="$peak" 'BEGIN { printf "%.1f", m / 1024 }')" \
    "$(glow_mem_total_gb)" "$(nproc)" "$rc" | tee -a "$GLOW_MON_FILE"
  rm -f "$GLOW_MON_PEAK" "$GLOW_MON_PEAK.tmp"
}
