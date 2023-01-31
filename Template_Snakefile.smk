# Template main snakefile for modular SuRE pipeline development.
import os
from git import Repo



# Establish snakefile paths.
SNAKEFILE = workflow.snakefile
REPO_DIR = os.path.dirname(os.path.abspath(SNAKEFILE))

#################################################################
# check status of code repository in git
# 1. snakemake can be run with option '--config noGitCheck=1' 
#    in this case the check is skipped
# 2. otherwise check git status and report commit hash or abort
if "noGitCheck" in config:
  print ("#######################################################")
  print ("  RUNNING WITHOUT CHECKING STATUS OF CODE REPOSITORY!  ")
  print ("#######################################################")
else:
  # get the Repo object corresponding to the pipeline git repo
  REPO = Repo(REPO_DIR)
  # assert it is not a bare repo (from examples)
  assert not REPO.bare

  # check if REPO is dirty (ie whether it contains modifications not yet commited)
  # if so, abort with message
  if REPO.is_dirty(untracked_files=True):
      print("\nYour snakefile is modified relative to the git repos", file=sys.stderr)
      print("Either undo your changes, or commit them", file=sys.stderr)
      print("Aborting\n\n", file=sys.stderr)
      sys.exit("")

  # Get and print the commit hash
  GIT_COMMIT_HASH = REPO.head.commit
  # print hash to stdout. This assumes output of the command snakemake is piped to a log file
  print ("git commit hash: ", GIT_COMMIT_HASH, "\n\n")
#################################################################
  

###########################################################################
###########################################################################
def clone_snakefile_module(MODULE, COMMIT, PATH_TO_SMK, VERBOSE=True):
# function to clone a module repository from github unless it is already present (and in correct status)

# The repo to clone may be already present, in which case we don't need to clone it again
# perform the following checks:
#   - check dir exists
#   - check dir is a git dir
#   - check repo is not a bare repo
#   - check repo is clean (no changes to code)
#   - check repo is in correct commit (by checking out the commit)
# If any of these checks fail the function either aborts the script or continues to clone the repo

  # define dir path where to clone the repo to
  REPOS="annogen/Snakefile_Modules"
  CLONE_TO=os.path.abspath(os.path.join(
      "modules",
      "Snakefile_Modules_{mod}_{com}".format(mod=MODULE, com=COMMIT)))
  PATH_TO_SMK=os.path.join(
      CLONE_TO,
      "Modules",
      PATH_TO_SMK)
  # Check dir to clone-to is already present locally
  if os.path.isdir(CLONE_TO):
    VERBOSE and print("repos dir ({}) exists".format(CLONE_TO))
    # dir is present; is it a git dir?
    try:
      _ = git.Repo(CLONE_TO).git_dir
      VERBOSE and print("repos dir is git dir")
      mod_repo = Repo(CLONE_TO)
      assert not mod_repo.bare
      VERBOSE and print("repos is not bare")
      # dir is a git repo; check repo is clean
      if not mod_repo.is_dirty(untracked_files=True):
        VERBOSE and print("repo not dirty")
        # repo is clean; check repo is at correct commit by checking it out
        try:
          mod_repo.git.checkout(COMMIT)
        except Exception as e:
          print("\ndirectory with cloned repo exists but can't checkout specified commit ({com})\n".format(com=COMMIT))
          sys.exit(e)
        # clean repo already exists locally, in correct commit; return
        VERBOSE and print("repo checked out clean")
        return(PATH_TO_SMK)
      else:
        print("\ndirectory with cloned repo exists but is not clean\n")
        sys.exit(1)
    except git.exc.InvalidGitRepositoryError:
      if len(os.listdir(CLONE_TO)) != 0:
        print("\ndirectory ({}) for cloning the module repo exists but is not a git repo, and not empty\n".format(CLONE_TO))
        sys.exit(1)
  else:
    VERBOSE and print("CLONE_TO ({dir}) does not exist; will clone from github".format(dir=CLONE_TO))
    
  # above checks allow control to reach to this point, meaning the repo needs to be cloned from github
  # get github token from the environment 
  TOKEN=os.getenv('GH_TOKEN', None)
  if TOKEN is None:
    print("\nno token for github found in environment; aborting\n")
    sys.exit(1)
  # the following URL allows to clone using a token
  URL='https://{}:x-oauth-basic@github.com/{}'.format(TOKEN, REPOS)
  try:
    repo=Repo.clone_from(URL, CLONE_TO)
  except Exception as e:
    print("\nfailed to clone module repo from: {url}; aborting\n".format(url=URL))
    sys.exit(e)
  # checkout the intended commit
  try:
    repo.git.checkout(COMMIT)
  except Exception as e:
    print("\nfailed to checkout out commit ({com}) from repo ({repo})".format(com=COMMIT,repo=REPOS))
    print("removing directory with cloned repo")
    # if checking out the commit failed we need to clean up the cloned dir
    try: 
      shutil.rmtree(CLONE_TO)
    except Exception as e:
      sys.exit(e)
    print("aborting\n")
    sys.exit(e)
  # cloning is succesful; return
  return(PATH_TO_SMK)
###########################################################################
###########################################################################

## EXAMPLES FOR INCLUDING MODULES: ##
##  ## MODULE Counts2coverage
##  MODULE="Counts2coverage"
##  COMMIT="896cc06"
##  PATH_TO_SMK = clone_snakefile_module(MODULE, COMMIT, "Counts2coverage/rules/counts2coverage.smk")
##  #PATH_TO_SMK="/data/home/ludo/repos/annogen/Snakefile_Modules/Modules/Counts2coverage/rules/counts2coverage.smk"
##  print(PATH_TO_SMK)
##  module counts2coverage:
##   snakefile: PATH_TO_SMK
##   config: config
##  
##  use rule * from counts2coverage as c2c_*
##  all_targets.extend(counts2coverage.counts2coverage_targets(''))
##  all_targets.extend(counts2coverage.counts2normScores_targets(''))
##  
##  
##  ## MODULE SuREcoverage2peaks
##  if "cDNA" in config and config['cDNA']['SAMPLES'] is not None:
##    MODULE="SuREcov2peaks"
##    COMMIT="896cc06"
##    PATH_TO_SMK = clone_snakefile_module(MODULE, COMMIT, "SuREcov2peaks/rules/SuREcoverage2peaks.smk")
##    #PATH_TO_SMK="/data/home/ludo/repos/annogen/Snakefile_Modules/Modules/SuREcov2peaks/rules/SuREcoverage2peaks.smk"
##    print(PATH_TO_SMK)
##    module coverage2peaks:
##      snakefile: PATH_TO_SMK
##      config: config
##    
##    use rule * from coverage2peaks as c2p_*
##    all_targets.extend(coverage2peaks.SuREcoverage2peaks_strandSpecific_targets())




#Load in the required parameters which will be strictly structured across all projects.
OUTDIR = config["OUTDIR"]
LIBRARY = config["LIBRARY"]
iPCR_SAMPLES = config["iPCR"]["SAMPLES"]
cDNA_SAMPLES = config["cDNA"]["SAMPLES"]
plDNA_SAMPLES = config["plDNA"]["SAMPLES"]


# This function returns the final targets of the pipeline based on the input section of the config file.
# These are the targets for the standard SuRE pipeline. Feel free to make adjustments based on your needs!
def get_all_targets():
    if iPCR_SAMPLES == None:
        ipcr_targets = []
    else:
        ipcr_targets = expand(os.path.join(OUTDIR, "SuRE-counts_{CHROM}.txt.gz"), CHROM=config["CHR_TARGET"])
    if cDNA_SAMPLES == None:
        cdna_targets = []
    else:
        cdna_targets = expand(os.path.join(OUTDIR, config["cDNA"]["OUTDIR"], "{SAMPLE}", "{SAMPLE}_trimmed_table.txt.gz"), SAMPLE=cDNA_SAMPLES.keys())
    if plDNA_SAMPLES == None:
        pldna_targets = []
    else:
        pldna_targets = expand(os.path.join(OUTDIR, config["plDNA"]["OUTDIR"], "{SAMPLE}", "{SAMPLE}_trimmed_table.txt.gz"), s=plDNA_SAMPLES.keys())
    return ipcr_targets + cdna_targets + pldna_targets

# Include all the snakefile modules that are required to produce the final targets.
# Repository structure ensures all modules are in the rules sub-directory relative to the snakefile.
include: os.path.join(REPO_DIR, "rules","iPCR_Trimming_Cutadapt.smk")


# Initialize the all rule to execute the included modules.
rule all:
    input:
        get_all_targets()
