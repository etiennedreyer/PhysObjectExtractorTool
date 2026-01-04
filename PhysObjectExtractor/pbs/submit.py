import os
import argparse

top_dir = "/srv01/agrp/dmitrykl/projects/cmssw/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor"
# top_out_dir = "/storage/agrp/dmitrykl/fastsim/f_delphes/data/full_event/raw_impact_all/"
top_out_dir = "/storage/agrp/dmitrykl/fastsim/f_delphes/data/full_event/raw_analysis/"
run_script = top_dir + "/pbs/run.sh"

# redo_list = [116, 149, 150, 21, 26, 47, 55, 94, 95]
redo_list = []

def submit_jobs(input, output=None, tag='tmp', limit="walltime=7:59:00,mem=12g,io=20", dry_run=False, num_events=-1, num_files=1000000):

    files = []
    if input.endswith(".root"):
        files.append(input)
    elif input.endswith(".txt"):
        with open(input) as f:
            files = f.readlines()
    else:
        raise ValueError("Input file must be a .root file or a .txt file containing a list of .root files")

    for jobID, file in enumerate(files):
        if jobID >= num_files:
            break
        jobID += 1
        file = file.strip()
        if file.startswith("#") or '.root' not in file:
            continue

        file = file.strip()
        if len(files)==1 and (output is not None):
            outfile = output
            out_dir = "/".join(outfile.split("/")[:-1])
        else:
            latter_part = file.split("/eos/opendata/cms/")[-1]
            out_dir = "/".join(latter_part.split("/")[:-1])
            out_dir = top_out_dir + '/' + out_dir
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            outfile = out_dir + "/" + file.split("/")[-1].replace(".root", "_poet.root")

        if len(redo_list)>0:
            if jobID not in redo_list:
                continue
        else:
            os.system(f"rm -rf {top_dir}/pbs/{tag}/out_{jobID}.log")
            os.system(f"rm -rf {top_dir}/pbs/{tag}/err_{jobID}.log")
            os.system(f"rm -rf {outfile}")

        os.makedirs(f"{top_dir}/pbs/{tag}", exist_ok=True)
        out_log = f"{top_dir}/pbs/{tag}/out_{jobID}.log"
        err_log = f"{top_dir}/pbs/{tag}/err_{jobID}.log"

        ### write job script in the output directory
        run_cms_script = f"{top_dir}/pbs/{tag}/run_cms_{jobID}.sh"
        with open(run_cms_script, "w") as f:
            # f.write(f"pwd\ncmsrel CMSSW_5_3_32\ncd CMSSW_5_3_32/src\ncmsenv\ncd PhysObjectExtractorTool/PhysObjectExtractor/\npwd\n")
            # f.write(f"cmsRun python/poet_cfg_genParticles.py {file} {outfile} {num_events}")
            # f.write(f'printf "blabla\\n" | CMS_PATH=~/projects/cmssw/ cmsRun python/poet_cfg_genParticles.py {file} {outfile} {num_events}')
            f.write(f'printf "blabla\\n" | CMS_PATH=~/projects/cmssw/ cmsRun python/poet_cfg_genParticles.py {file} {outfile} {num_events}')
        job_name = f"{tag}_{jobID}"
        job = f"qsub -N {job_name} -q N -l {limit} -o {out_log} -e {err_log} -v tag={tag},jobID={jobID} {run_script}"
        print(job)
        if not dry_run:
            os.system(job)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Submit poet jobs to the PBS batch scheduler")
    parser.add_argument("-i", "--input", required=True, help="Input file or list of input files")
    parser.add_argument("-o", "--output", help="Output file", default=None)
    parser.add_argument("-t", "--tag", help="Tag for temp job files", default="tmp")
    parser.add_argument("-l", "--limit", help="Resource limit", default="walltime=4:59:00,mem=2g,io=2")
    parser.add_argument("-d", "--dry_run", help="Dry run", action="store_true")
    parser.add_argument("-n", "--num_events", help="Number of events to process", default=1000000)
    parser.add_argument("-nf", "--num_files", help="Number of files to process", default=1000000, type=int)
    args = parser.parse_args()

    submit_jobs(args.input, args.output, tag=args.tag, limit=args.limit, dry_run=args.dry_run, num_events=args.num_events, num_files=args.num_files)