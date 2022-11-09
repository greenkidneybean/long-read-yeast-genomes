#!/usr/bin/env bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --time=3:59:59
#SBATCH --gres=lscratch:100
#SBATCH --partition=quick,norm
#SBATCH --output %j.slurm.out
#SBATCH --error %j.slurm.out

## PARSE ARGUMENTS
# Convenience functions.
usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }

# One loop, nothing more.
if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in

            # Your options go here.
            -h|--help) HELP='true';;
            -f|--flag) flag='true';;
            --strain1) assert_argument "$1" "$opt"; STRAIN1="$1"; shift;;
            --strain2) assert_argument "$1" "$opt"; STRAIN2="$1"; shift;;
            --chr) assert_argument "$1" "$opt"; CHR="$1"; shift;;
            --reference) assert_argument "$1" "$opt"; REF="$1"; shift;;
      
            # Arguments processing. You may remove any unneeded line after the 1st.
            -|''|[!-]*) set -- "$@" "$opt";;                                          # positional argument, rotate to the end
            --*=*)      set -- "${opt%%=*}" "${opt#*=}" "$@";;                        # convert '--name=arg' to '--name' 'arg'
            -[!-]?*)    set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@";;       # convert '-abc' to '-a' '-b' '-c'
            --)         while [ "$1" != "$EOL" ]; do set -- "$@" "$1"; shift; done;;  # process remaining arguments as positional
            -*)         usage_error "unknown option: '$opt'";;                        # catch misspelled options
            *)          usage_error "this should NEVER happen ($opt)";;               # sanity test for previous patterns
    
        esac
    done
    shift  # $EOL
fi

if [ "${HELP}" == 'true' ]; then

cat << EndOfHelp
    --help | -h         show this message
    --strain1           first strain in cross (required)
    --strain2           second strain in cross (required)
    --chr               chromosome being mapped (required)
    --reference         fasta file for mapping (required)
EndOfHelp
exit 0
fi


argexit='false'
if [ "${CHR}" == '' ]; then echo "ERROR: --chr is required"; argexit='true'
else echo "--chromosome is ${CHR}"; fi

if [ "${STRAIN1}" == '' ]; then echo "ERROR: --strain1 is required"; argexit='true'
else echo "--strain1 is ${STRAIN1}"; fi

if [ "${STRAIN2}" == '' ]; then echo "ERROR: --strain2 is required"; argexit='true'
else echo "--strain2 is ${STRAIN2}"; fi

if [ "${REF}" == '' ]; then echo "ERROR: --reference is required"; argexit='true'
else echo "--reference is ${REF}"; fi

if [ ${argexit} == 'true' ]; then
    echo 'exiting due to missing arguments'
    exit 1
fi

## END ARGUMENT PARSING

