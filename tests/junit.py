#!/usr/bin/env python3

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'junit-xml')))
from junit_xml import TestCase, TestSuite

def parse_testargs(file):
    allargs = []
    if os.path.splitext(file)[1] in ['.c', '.cpp']:
        for line in open(file).readlines():
            if line.startswith('//TESTARGS'):
                args = [line.split()[1:]]
                args.append([line.split()[0].strip('//TESTARGS(name=').strip(')')])
                allargs.append(args)
        return allargs
    elif os.path.splitext(file)[1] == '.usr':
        for line in open(file).readlines():
            if line.startswith('C_TESTARGS'):
                args = [line.split()[1:]]
                args.append([line.split()[0].strip('C_TESTARGS(name=').strip(')')])
                allargs.append(args)
        return allargs
    raise RuntimeError('Unrecognized extension for file: {}'.format(file))

def get_source(test):
    if test.startswith('petsc-'):
        return os.path.join('examples', 'petsc', test[6:] + '.c')
    elif test.startswith('mfem-'):
        return os.path.join('examples', 'mfem', test[5:] + '.cpp')
    elif test.startswith('nek-'):
        return os.path.join('examples', 'nek', 'bps', test[4:] + '.usr')
    elif test.startswith('fluids-'):
        return os.path.join('examples', 'fluids', test[7:] + '.c')
    elif test.startswith('solids-'):
        return os.path.join('examples', 'solids', test[7:] + '.c')
    elif test.startswith('ex'):
        return os.path.join('examples', 'ceed', test + '.c')

def get_testargs(test):
    source = get_source(test)
    if source is None:
        return [[['{ceed_resource}'], ['']]]
    return parse_testargs(source)

def check_required_failure(case, stderr, required):
    if required in stderr:
        case.status = 'fails with required: {}'.format(required)
    else:
        case.add_failure_info('required: {}'.format(required))

def contains_any(resource, substrings):
    return any((sub in resource for sub in substrings))

def skip_rule(test, resource):
    return any((
        test.startswith('fluids-') and contains_any(resource, ['occa', 'gpu']) and not contains_any(resource, ['/gpu/cuda/gen']),
        test.startswith('solids-') and contains_any(resource, ['occa']),
        test.startswith('nek') and contains_any(resource, ['occa']),
        test.startswith('t507') and contains_any(resource, ['occa']),
        test.startswith('t318') and contains_any(resource, ['magma']),
        test.startswith('t506') and contains_any(resource, ['magma']),
        ))
        
def run(test, backends):
    import subprocess
    import time
    import difflib
    allargs = get_testargs(test)

    testcases = []
    for args, name in allargs:
        for ceed_resource in backends:
            rargs = [os.path.join('build', test)] + args.copy()
            rargs[rargs.index('{ceed_resource}')] = ceed_resource

            if skip_rule(test, ceed_resource):
                case = TestCase('{} {}'.format(test, ceed_resource),
                                elapsed_sec=0,
                                timestamp=time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime()),
                                stdout='',
                                stderr='')
                case.add_skipped_info('Pre-run skip rule')
            else:
                start = time.time()
                proc = subprocess.run(rargs,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
                proc.stdout = proc.stdout.decode('utf-8')
                proc.stderr = proc.stderr.decode('utf-8')

                case = TestCase('{} {} {}'.format(test, *name, ceed_resource),
                                elapsed_sec=time.time()-start,
                                timestamp=time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(start)),
                                stdout=proc.stdout,
                                stderr=proc.stderr)
                ref_stdout = os.path.join('tests/output', test + '.out')

            if not case.is_skipped() and proc.stderr:
                if 'OCCA backend failed to use' in proc.stderr:
                    case.add_skipped_info('occa mode not supported {} {}'.format(test, ceed_resource))
                elif 'Backend does not implement' in proc.stderr:
                    case.add_skipped_info('not implemented {} {}'.format(test, ceed_resource))
                elif 'Can only provide to HOST memory' in proc.stderr:
                    case.add_skipped_info('device memory not supported {} {}'.format(test, ceed_resource))

            if not case.is_skipped():
                if test[:4] in 't110 t111 t112 t113 t114'.split():
                    check_required_failure(case, proc.stderr, 'Cannot grant CeedVector array access')
                if test[:4] in 't115'.split():
                    check_required_failure(case, proc.stderr, 'Cannot grant CeedVector read-only array access, the access lock is already in use')
                if test[:4] in 't116'.split():
                    check_required_failure(case, proc.stderr, 'Cannot destroy CeedVector, the writable access lock is in use')
                if test[:4] in 't117'.split():
                    check_required_failure(case, proc.stderr, 'Cannot restore CeedVector array access, access was not granted')
                if test[:4] in 't118'.split():
                    check_required_failure(case, proc.stderr, 'Cannot sync CeedVector, the access lock is already in use')
                if test[:4] in 't215'.split():
                    check_required_failure(case, proc.stderr, 'Cannot destroy CeedElemRestriction, a process has read access to the offset data')
                if test[:4] in 't303'.split():
                    check_required_failure(case, proc.stderr, 'Length of input/output vectors incompatible with basis dimensions')

            if not case.is_skipped() and not case.status:
                if proc.stderr:
                    case.add_failure_info('stderr', proc.stderr)
                elif proc.returncode != 0:
                    case.add_error_info('returncode = {}'.format(proc.returncode))
                elif os.path.isfile(ref_stdout):
                    with open(ref_stdout) as ref:
                        diff = list(difflib.unified_diff(ref.readlines(),
                                                         proc.stdout.splitlines(keepends=True),
                                                         fromfile=ref_stdout,
                                                         tofile='New'))
                    if diff:
                        case.add_failure_info('stdout', output=''.join(diff))
                elif proc.stdout and test[:4] not in 't003':
                    case.add_failure_info('stdout', output=proc.stdout)
            testcases.append(case)
    return TestSuite(test, testcases)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Test runner with JUnit output')
    parser.add_argument('--output', help='Output file to write test', default=None)
    parser.add_argument('--gather', help='Gather all *.junit files into XML', action='store_true')
    parser.add_argument('test', help='Test executable', nargs='?')
    args = parser.parse_args()

    if args.gather:
        gather()
    else:
        backends = os.environ['BACKENDS'].split()

        result = run(args.test, backends)
        output = (os.path.join('build', args.test + '.junit')
                  if args.output is None
                  else args.output)
        with open(output, 'w') as fd:
            TestSuite.to_file(fd, [result])
        for t in result.test_cases:
            failures = len([c for c in result.test_cases if c.is_failure()])
            errors = len([c for c in result.test_cases if c.is_error()])
            if failures + errors > 0:
                sys.exit(1)
