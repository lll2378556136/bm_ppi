#!/opt/conda/bin/python
# Patched copy of DPAM v2.0 step8_iterative_dali.py (upstream:
# CongLabCode/DPAM docker/scripts/step8_iterative_dali.py, md5
# febb71be2624cbec4d19549cd752ad0e). scripts/run_dpam.slurm bind-mounts
# this file over /opt/DPAM/scripts/step8_iterative_dali.py in sifs/dpam.sif.
#
# Why patched: upstream workers call sys.exit(1) when dali.pl fails. A worker
# killed by sys.exit never resolves its ApplyResult, so the parent blocked
# forever in process.get() (observed on protein 101743456: 8336 SIGSEGV
# crashes of DSSP/puutos, deadlocked job idle >11h). Changes vs upstream:
#   1. run_dali reports dali.pl failure as a return value and cleans its /tmp
#      scratch; unexpected exceptions are converted to task failures too.
#   2. Each task starts from a fresh /tmp work dir (stale-scratch guard).
#   3. The parent polls results with a failure-fraction circuit breaker and a
#      per-protein wall-clock cap instead of a blocking get(). An aborted
#      protein keeps the hits that finished; if none did, it gets the standard
#      .done marker, which run_step8.py treats as a terminal zero-hit state.
import os, sys
import time, subprocess, traceback
from multiprocessing import Pool

wdir = os.getcwd()

# Circuit-breaker tunables (env-overridable): abort a protein once at least
# FAIL_MIN tasks finished and FAIL_FRACTION of them failed; WALL_TIMEOUT_S is
# the per-protein backstop against tasks that never finish (e.g. hung dali.pl).
FAIL_MIN = int(os.environ.get('DPAM_STEP8_FAIL_MIN', '64'))
FAIL_FRACTION = float(os.environ.get('DPAM_STEP8_FAIL_FRAC', '0.95'))
WALL_TIMEOUT_S = int(os.environ.get('DPAM_STEP8_WALL_TIMEOUT', str(6 * 3600)))
POLL_INTERVAL_S = int(os.environ.get('DPAM_STEP8_POLL_INTERVAL', '15'))

def get_domain_range(resids):
    segs = []
    resids.sort()
    cutoff1 = 5
    cutoff2 = len(resids) * 0.05
    cutoff = max(cutoff1, cutoff2)
    for resid in resids:
        if not segs:
            segs.append([resid])
        else:
            if resid > segs[-1][-1] + cutoff:
                segs.append([resid])
            else:
                segs[-1].append(resid)
    seg_string = []
    for seg in segs:
        start = str(seg[0])
        end = str(seg[-1])
        seg_string.append(start + '-' + end)
    return ','.join(seg_string)

def cleanup_scratch(prot, edomain):
    # Remove one task's /tmp scratch; <prot>_<edomain>_hits is deliberately
    # kept so hits salvaged from earlier iterations still reach the count.
    # chdir out first so the next os.system does not start with a deleted cwd.
    try:
        os.chdir('/tmp/')
    except Exception:
        pass
    os.system('rm -rf /tmp/tmp_' + prot + '_' + edomain)
    os.system('rm -f /tmp/' + prot + '_' + edomain + '.pdb')

def run_dali(edomain):
    # Return 0 on success, 1 on dali.pl failure or unexpected exception.
    try:
        os.system('cp ' + wdir + '/step2/' + dataset + '/' + prot + '.pdb /tmp/' + prot + '_' + edomain + '.pdb')
        alicount = 0
        # Start from a fresh work dir: scratch left by a crashed or killed
        # attempt used to leak into the next run of the same pair.
        os.system('rm -rf /tmp/tmp_' + prot + '_' + edomain)
        os.system('mkdir /tmp/tmp_' + prot + '_' + edomain)
        os.chdir('/tmp/tmp_' + prot + '_' + edomain)
        while 1:
            rcode=subprocess.run('dali.pl --pdbfile1 /tmp/' + prot + '_' + edomain + '.pdb --pdbfile2 /mnt/databases/ECOD_pdbs/' + edomain + '.pdb --dat1 ./ -dat2 ./ --outfmt summary,alignments,transrot > log',shell=True).returncode
            if str(rcode) != '0':
                # Upstream sys.exit(1)'d here, killing the worker without
                # resolving its ApplyResult and deadlocking the parent.
                cleanup_scratch(prot, edomain)
                return 1
            fp = os.popen('ls -1 mol*.txt')
            info = fp.readlines()
            fp.close()
            filenames = []
            for line in info:
                filenames.append(line[:-1])

            if filenames:
                fp = open('/tmp/' + prot + '_' + edomain + '.pdb', 'r')
                Qresids_set = set([])
                for line in fp:
                    resid = int(line[22:26])
                    Qresids_set.add(resid)
                fp.close()
                Qresids = list(Qresids_set)

                info = []
                for filename in filenames:
                    fp = open(filename, "r")
                    lines = fp.readlines()
                    fp.close()
                    for line in lines:
                        info.append(line)

                qali = ''
                sali = ''
                rotdata = []
                getit = 1
                zscore = 0
                for line in info:
                    words = line.split()
                    if len(words) >= 2 and getit:
                        if words[0] == 'Query':
                            qali += words[1]
                        elif words[0] == 'Sbjct':
                            sali += words[1]
                        elif words[0] == 'No' and words[1] == '1:':
                            for word in words:
                                if '=' in word:
                                    subwords = word.split('=')
                                    if subwords[0] == 'Z-score':
                                        zinfo = subwords[1].split('.')
                                        zscore = float(zinfo[0] + '.' + zinfo[1])
                        elif len(words) >= 8 and words[0] == "-matrix":
                            rotdata.append(line)
                        elif words[0] == 'No' and words[1] == '2:':
                            getit = 0

                qinds = []
                sinds = []
                length = len(qali)
                qposi = 0
                sposi = 0
                match = 0
                for i in range(length):
                    if qali[i] != '-':
                        qposi += 1
                    if sali[i] != '-':
                        sposi += 1
                    if qali[i] != '-' and sali[i] != '-':
                        if qali[i].isupper() and sali[i].isupper():
                            match += 1
                            qinds.append(qposi)
                            sinds.append(sposi)
                qlen = qposi
                slen = sposi

                if match >= 20:
                    alicount += 1
                    rp = open('/tmp/' + prot + '_' + edomain + '_hits', 'a')
                    rp.write('>' + edomain + '_' + str(alicount) + '\t' + str(zscore) + '\t' + str(match) + '\t' + str(qlen) + '\t' + str(slen) + '\n')
                    rotation = ''
                    translation = 'translation'
                    for rotline in rotdata:
                        rotwords = rotline[:-2].split()
                        rotation += 'rotation\t' + '\t'.join(rotwords[4:7]) + '\n'
                        translation += '\t' + rotwords[7]
                    rp.write(rotation)
                    rp.write(translation + '\n')
                    for i in range(len(qinds)):
                        qind = qinds[i] - 1
                        sind = sinds[i]
                        rp.write(str(Qresids[qind]) + '\t' + str(sind) + '\n')
                    rp.close()

                    raw_qresids = []
                    for qind in qinds:
                        raw_qresids.append(Qresids[qind - 1])
                    qrange = get_domain_range(raw_qresids)
                    qresids = set([])
                    qsegs = qrange.split(',')
                    for qseg in qsegs:
                        qedges = qseg.split('-')
                        qstart = int(qedges[0])
                        qend = int(qedges[1])
                        for qres in range(qstart, qend + 1):
                            qresids.add(qres)
                    remain_resids = Qresids_set.difference(qresids)

                    if len(remain_resids) >= 20:
                        rp = open('/tmp/' + prot + '_' + edomain + '.pdbnew', 'w')
                        fp = open('/tmp/' + prot + '_' + edomain + '.pdb', 'r')
                        for line in fp:
                            resid = int(line[22:26])
                            if resid in remain_resids:
                                rp.write(line)
                        fp.close()
                        rp.close()
                        os.system('mv /tmp/' + prot + '_' + edomain + '.pdbnew /tmp/' + prot + '_' + edomain + '.pdb')
                        if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                            pass
                            os.system('rm *')
                    else:
                        if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                            pass
                            os.system('rm *')
                        break
                else:
                    if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                        pass
                        os.system('rm *')
                    break
            else:
                if os.getcwd() == '/tmp/tmp_' + prot + '_' + edomain:
                    pass
                    os.system('rm *')
                break
        os.chdir('/tmp/')
        os.system('rmdir /tmp/tmp_' + prot + '_' + edomain)
        os.system('rm /tmp/' + prot + '_' + edomain + '.pdb')
        return 0
    except Exception:
        # An escaping exception used to propagate through ApplyResult.get()
        # and crash the whole per-protein script; report it as a task failure.
        traceback.print_exc()
        cleanup_scratch(prot, edomain)
        return 1


dataset = sys.argv[1]
prot = sys.argv[2]
ncore = int(sys.argv[3])

if os.path.exists(wdir + '/step7/' + dataset + '/' + prot + '_hits'):
    fp = open(wdir + '/step7/' + dataset + '/' + prot + '_hits','r')
    edomains = []
    for line in fp:
        words = line.split()
        edomains.append(words[0])
    fp.close()

    inputs = []
    for edomain in edomains:
        inputs.append([edomain])
    pool = Pool(processes = ncore)
    results = []
    for item in inputs:
        process = pool.apply_async(run_dali, item)
        results.append(process)

    # Upstream blocked here with `for process in results: process.get()`.
    # Poll instead so systematic dali.pl failures are counted and cut short.
    aborted = False
    n_done = 0
    n_fail = 0
    t0 = time.time()
    while 1:
        n_done = 0
        n_fail = 0
        for process in results:
            if process.ready():
                n_done += 1
                if process.get():
                    n_fail += 1
        if n_done == len(results):
            break
        if n_done >= FAIL_MIN and n_fail >= FAIL_FRACTION * n_done:
            aborted = True
            break
        if time.time() - t0 > WALL_TIMEOUT_S:
            aborted = True
            break
        time.sleep(POLL_INTERVAL_S)

    if aborted:
        sys.stderr.write('step8_iterative_dali: aborted %s after %d/%d tasks finished (%d failed, %d seconds elapsed)\n' % (prot, n_done, len(results), n_fail, int(time.time() - t0)))
        sys.stderr.flush()
        pool.terminate()
        pool.join()
        # Scratch left by killed workers; *_hits files are kept for the count.
        os.system('rm -rf /tmp/tmp_' + prot + '_*')
        os.system('rm -f /tmp/' + prot + '_*.pdb')
    else:
        pool.close()
        pool.join()

    fp = os.popen('ls -1 /tmp/' + prot + '_*_hits')
    get_hit_count = 0
    for line in fp:
        get_hit_count += 1
    fp.close()
    if get_hit_count:
        os.system('cat /tmp/' + prot + '_*_hits > ' + wdir + '/step8/' + dataset + '/' + prot + '_hits')
        os.system('rm /tmp/' + prot + '_*_hits')
    else:
        os.system('echo \'done\' > ' + wdir + '/step8/' + dataset + '/' + prot + '.done')
else:
    os.system('echo \'done\' > ' + wdir + '/step8/' + dataset + '/' + prot + '.done')
