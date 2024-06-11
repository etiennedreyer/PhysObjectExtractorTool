import os
import sys
import argparse

top_dir = "/storage/agrp/dreyet/CMS_OpenData/FullEvent/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor"
run_script = top_dir + "/pbs/run.sh"

def submit_jobs(input, output=None, limit="walltime=7:59:00,mem=12g", dry_run=False, num_events=-1):

    files = []
    if input.endswith(".root"):
        files.append(input)
    elif input.endswith(".txt"):
        with open(input) as f:
            files = f.readlines()
    else:
        raise ValueError("Input file must be a .root file or a .txt file containing a list of .root files")

    for file in files:
        if len(files)==1 and (output is not None):
            outfile = output
        else:
            latter_part = file.split("/eos/opendata/cms/")[-1]
            out_dir = latter_part.split("/")[:-1].join("/")
            out_dir = top_dir + "/output/" + out_dir
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            outfile = out_dir + "/" + file.split("/")[-1].replace(".root", "_poet.root")

        out_log = outfile.replace(".root", "_poet_out.log")
        err_log = outfile.replace(".root", "_poet_err.log")

        job = f"qsub -q N -l {limit} -o {out_log} -e {err_log} -v infile={file},outfile={outfile},nevents={num_events} {run_script}"
        print(job)
        if not dry_run:
            os.system(job)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Submit poet jobs to the PBS batch scheduler")
    parser.add_argument("-i", "--input", required=True, help="Input file or list of input files")
    parser.add_argument("-o", "--output", help="Output file", default=None)
    parser.add_argument("-l", "--limit", help="Resource limit", default="walltime=7:59:00,mem=12g")
    parser.add_argument("-d", "--dry_run", help="Dry run", action="store_true")
    parser.add_argument("-n", "--num_events", help="Number of events to process", default=-1)
    args = parser.parse_args()

    submit_jobs(args.input, args.output, limit=args.limit, dry_run=args.dry_run, num_events=args.num_events)