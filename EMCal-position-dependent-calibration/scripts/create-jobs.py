#!/usr/bin/env python3
import numpy as np
import subprocess
import argparse
import os
import shutil

parser = argparse.ArgumentParser()
subparser = parser.add_subparsers(dest='command')

create = subparser.add_parser('create')
status = subparser.add_parser('status')
hadd   = subparser.add_parser('hadd')

create.add_argument('--executable', type=str, help='Job script to execute.', required=True)
create.add_argument('--macros', type=str, help='Directory of input macros. Directory containing Fun4All_G4_sPHENIX.C and G4Setup_sPHENIX.C.',required=True)
create.add_argument('-n', '--events', type=int, default=1, help='Number of events to generate. Default: 1.')
create.add_argument('-d', '--output', type=str, default='test', help='Output Directory. Default: Current Directory.')
create.add_argument('-m', '--jobs-per-submission', type=int, default=20000, help='Maximum number of jobs per condor submission. Default: 20000.')
create.add_argument('-j', '--events-per-job', type=int, default=50, help='Number of events to generate per job. Default: 50.')
create.add_argument('--memory', type=int, default=10, help='Memory (units of GB) to request per condor submission. Default: 10 GB.')

status.add_argument('-d','--condor-dir', type=str, help='Condor submission directory.', required=True)

hadd.add_argument('-i','--job-dir-list', type=str, help='List of directories containing condor jobs to be merged.', required=True)
hadd.add_argument('-o','--output', type=str, default='test.root', help='Output root file. Default: test.root.')
hadd.add_argument('-n','--jobs-per-hadd', type=int, default=5000, help='Number of jobs to merge per hadd call. Default: 5000.')
hadd.add_argument('-j','--jobs-open', type=int, default=50, help='Number of jobs to load at once. Default: 50.')

args = parser.parse_args()

def create_jobs():
    events              = args.events
    jobs_per_submission = args.jobs_per_submission
    output_dir          = os.path.abspath(args.output)
    executable          = os.path.abspath(args.executable)
    events_per_job      = args.events_per_job
    memory              = args.memory
    macros_dir          = os.path.abspath(args.macros)
    jobs                = events//events_per_job
    submissions         = int(np.ceil(jobs/jobs_per_submission))

    print(f'Events: {events}')
    print(f'Events per job: {events_per_job}')
    print(f'Jobs: {jobs}')
    print(f'Maximum jobs per condor submission: {jobs_per_submission}')
    print(f'Submissions: {submissions}')
    print(f'Requested memory per job: {memory}GB')
    print(f'Output Directory: {output_dir}')
    print(f'Macros Directory: {macros_dir}')
    print(f'Executable: {executable}')
    print('-----------------------------------')

    os.makedirs(output_dir,exist_ok=True)
    # Generate condor submission file
    condor_file = f'{output_dir}/genFun4All.sub'
    with open(condor_file, mode="w") as file:
        file.write(f'executable             = bin/{os.path.basename(executable)}\n')
        file.write(f'arguments              = $(myPid) $(initialSeed) {events_per_job}\n')
        file.write('log                     = log/job-$(myPid).log\n')
        file.write('output                  = stdout/job-$(myPid).out\n')
        file.write('error                   = error/job-$(myPid).err\n')
        file.write(f'request_memory         = {memory}GB\n')
        file.write('should_transfer_files   = YES\n')
        file.write('when_to_transfer_output = ON_EXIT\n')
        file.write('transfer_input_files    = src/Fun4All_G4_sPHENIX.C, src/G4Setup_sPHENIX.C\n')
        file.write('transfer_output_files   = G4sPHENIX_g4cemc_eval-$(myPid).root\n')
        file.write('transfer_output_remaps  = "G4sPHENIX_g4cemc_eval-$(myPid).root = output/G4sPHENIX_g4cemc_eval-$(myPid).root"\n')
        file.write('queue myPid,initialSeed from seed.txt')

    for i in range(submissions):
        print(f'Submission: {i}')

        submit_dir = f'{output_dir}/submission-{i}'
        print(f'Submission Directory: {submit_dir}')

        os.makedirs(f'{submit_dir}/stdout',exist_ok=True)
        os.makedirs(f'{submit_dir}/error',exist_ok=True)
        os.makedirs(f'{submit_dir}/log',exist_ok=True)
        os.makedirs(f'{submit_dir}/output',exist_ok=True)
        os.makedirs(f'{submit_dir}/bin',exist_ok=True)
        os.makedirs(f'{submit_dir}/src',exist_ok=True)

        shutil.copy(condor_file, submit_dir)
        shutil.copy(executable, f'{submit_dir}/bin')
        shutil.copy(f'{macros_dir}/Fun4All_G4_sPHENIX.C', f'{submit_dir}/src')
        shutil.copy(f'{macros_dir}/G4Setup_sPHENIX.C', f'{submit_dir}/src')

        file_name = f'{submit_dir}/seed.txt'
        with open(file_name, mode="w") as file:
            for j in range(min(jobs,jobs_per_submission)):
                file.write(f'{j} {i*jobs_per_submission+100}\n')

        jobs -= min(jobs,jobs_per_submission)
        print(f'Written {file_name}')

def get_status():
    condor_dir = os.path.abspath(args.condor_dir)
    submit_dirs = next(os.walk(condor_dir))[1]
    print(f'Condor Directory: {condor_dir}')
    jobs_done_total = 0
    total = 0
    for submit_dir in submit_dirs:
        jobs_done = len(os.listdir(f'{condor_dir}/{submit_dir}/output'))
        jobs_total = len(os.listdir(f'{condor_dir}/{submit_dir}/log'))
        if(jobs_total != 0):
            print(f'Condor submission dir: {condor_dir}/{submit_dir}, done: {jobs_done}, {jobs_done/jobs_total*100:.2f} %')
            jobs_done_total += jobs_done
            total += jobs_total

    if(total != 0):
        print(f'Total jobs done: {jobs_done_total}, {jobs_done_total/total*100:.2f} %')

def hadd():
    job_dir_list  = os.path.abspath(args.job_dir_list)
    output        = os.path.abspath(args.output)
    jobs_per_hadd = args.jobs_per_hadd
    jobs_open     = args.jobs_open+1
    print(f'file list: {job_dir_list}')
    print(f'output: {output}')
    print(f'jobs per hadd: {jobs_per_hadd}')
    print(f'jobs open at once: {jobs_open-1}')

    jobs = []
    with open(job_dir_list) as f:
        for line in f:
            line = line.strip()
            jobs_l = os.listdir(line)
            print(f'dir: {line}, jobs: {len(jobs_l)}')
            jobs.extend([os.path.join(line,file) for file in jobs_l])

    total_jobs = len(jobs)
    hadd_calls = int(np.ceil(total_jobs/jobs_per_hadd))

    print(f'total jobs: {total_jobs}')
    print(f'hadd calls: {hadd_calls}')

    for i in range(hadd_calls):
        subprocess.run(['echo', '#######################'])
        subprocess.run(['echo', f'working on hadd: {i}'])
        command = ['hadd', '-a', '-n', str(jobs_open), output]
        i_start = jobs_per_hadd*i
        i_end = min(jobs_per_hadd*(i+1), total_jobs)
        subprocess.run(['echo', f'i_start: {i_start}, i_end: {i_end}'])
        command.extend(jobs[i_start:i_end])
        subprocess.run(command)
        subprocess.run(['echo', f'done with hadd: {i}'])
        subprocess.run(['echo', '#######################'])

if __name__ == '__main__':
    if(args.command == 'create'):
        create_jobs()
    elif(args.command == 'status'):
        get_status()
    elif(args.command == 'hadd'):
        hadd()
