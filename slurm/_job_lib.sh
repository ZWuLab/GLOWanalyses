# ============================================================================
# slurm/_job_lib.sh - shared helpers for the slurm/run-*.sh array templates
#
# SOURCED (not run) by every slurm/run-*.sh, right after its CWD guard, from the
# GLOWanalyses directory (`source slurm/_job_lib.sh`). Two helpers:
#   glow_require_config <config> <usage-hint>   fail fast, with the path + CWD,
#                                               when the run config is missing
#   glow_activate_r                             put the GLOW R environment on PATH
#
# R ENVIRONMENT INSIDE A JOB. A batch job inherits only the submit shell's
# environment variables; the R environment holding GFisher/GLOWr/GLOWpipeline
# must be activated here, inside the job. Two knobs, read as environment
# variables (export them before `sbatch`, or pass them on the sbatch line, e.g.
# `sbatch --export=ALL,GLOW_CONDA_ENV=GLOW ...`; edit the default below to make
# a choice permanent for your clone):
#   GLOW_CONDA_ENV   conda env (name or full path) to activate. Default: r_env.
#   GLOW_RSCRIPT     full path to an Rscript to use INSTEAD of conda, e.g. a
#                    `module load R` install, or .../envs/<env>/bin/Rscript.
#                    Its bin/ dir is prepended to PATH, so sibling command-line
#                    tools installed next to it (xsv, plink) resolve too.
#
# WHY nounset IS RELAXED AROUND `conda activate`: conda runs per-package
# (de)activation hooks (e.g. the gcc_linux-64 / gxx / gfortran compiler hooks)
# that reference variables which may be unset (SYS_SYSROOT,
# _CONDA_PYTHON_SYSCONFIGDATA_NAME_USED, ...). Under `set -u` such a reference
# aborts the job before R starts -- "<hook>.sh: line N: VAR: unbound variable" --
# for an env that activates fine interactively. nounset is switched off for the
# activation only and restored right after.

# glow_require_config <config> <usage-hint>
# The run config must be given (first positional arg, or $GLOW_CONFIG) and exist.
# Paths resolve from the GLOWanalyses directory (the job CWD), so the message
# prints the CWD -- a typo'd config path is the most common submit mistake.
glow_require_config() {
  local cfg="$1" usage="$2"
  if [ -z "${cfg}" ]; then
    echo "ERROR: pass the run config as the first arg (or set GLOW_CONFIG), e.g." >&2
    echo "  ${usage}" >&2
    return 1
  fi
  if [ ! -f "${cfg}" ]; then
    echo "ERROR: run config not found: ${cfg}  (CWD: ${PWD})" >&2
    echo "       Paths resolve from the GLOWanalyses directory: pass an absolute path," >&2
    echo "       or one relative to it (e.g. runs/<name>/config.R)." >&2
    return 1
  fi
}

# glow_activate_r
# Make `Rscript` (with the GLOW packages) available: GLOW_RSCRIPT if set, else
# `conda activate $GLOW_CONDA_ENV` (default r_env). Prints the Rscript used.
glow_activate_r() {
  if [ -n "${GLOW_RSCRIPT:-}" ]; then
    if [ ! -x "${GLOW_RSCRIPT}" ]; then
      echo "ERROR: GLOW_RSCRIPT is not an executable file: ${GLOW_RSCRIPT}" >&2
      return 1
    fi
    PATH="$(cd "$(dirname "${GLOW_RSCRIPT}")" && pwd):${PATH}"
    export PATH
  else
    local env="${GLOW_CONDA_ENV:-r_env}"
    if ! command -v conda >/dev/null 2>&1; then
      echo "ERROR: conda is not on PATH inside the job, so env '${env}' cannot be activated." >&2
      echo "       Submit from a shell where conda is initialised (sbatch exports its" >&2
      echo "       environment), or set GLOW_RSCRIPT=/path/to/Rscript to skip conda." >&2
      return 1
    fi
    local nounset=0
    case "$-" in *u*) nounset=1 ;; esac
    set +u      # conda's (de)activation hooks are not nounset-clean (see header)
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    if ! conda activate "${env}"; then
      if [ "${nounset}" = 1 ]; then set -u; fi
      echo "ERROR: 'conda activate ${env}' failed." >&2
      echo "       Set GLOW_CONDA_ENV to the conda env (name or path) holding the GLOW" >&2
      echo "       packages, or GLOW_RSCRIPT=/path/to/Rscript to skip conda." >&2
      return 1
    fi
    if [ "${nounset}" = 1 ]; then set -u; fi
  fi
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript is not on PATH after activating the R environment." >&2
    return 1
  fi
  echo "R environment: $(command -v Rscript)"
}
