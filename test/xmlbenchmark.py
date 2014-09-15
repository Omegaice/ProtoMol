import glob
import os
import subprocess
import optparse
import csv

# Arguments
parser = optparse.OptionParser(description='ProtoMol Test Suite')
parser.add_option('--verbose', '-v', action='store_true', dest='verbose', default=False, help='Verbose output')
parser.add_option('--parallel', '-p', action='store_true', dest='parallel', default=False, help='MPI Testing')

(options, args) = parser.parse_args()

# Empty Output Directory
outdata = glob.glob( os.path.join(os.getcwd(),'tests','output','*' ) )
map( os.remove, outdata )

# Find Executable
executable = []
if options.parallel:
    executable.append( "mpirun" )
    executable.append( "-np" )
    executable.append( "4" )
executable.append( os.path.join(os.getcwd(), 'ProtoMol') )

# Find Tests
tests = glob.glob(os.path.join(os.getcwd(), 'tests', '*-benchmark.conf'))
tests.sort()

stats_test = len(tests)

# Run Tests
times = []

for test in tests:
    name = os.path.splitext(os.path.basename(test))[0]

    # Run ProtoMol
    cmd = executable[:]
    cmd.append( test )

    print "Executing Benchmark: " + name
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()

    time = 0.0

    lines = stderr.split('\n')
    for line in lines:
        if line.find("Timing: ") != -1:
            data = line.split()
            time = float(data[2].replace('[s]', ''))

    print "Time: " + str(time) + "ms"
    print ""

    times.append( time )


# Write CSV
writer = csv.writer(open("performance.csv", "wb"))

def shortname( path ):
    return os.path.splitext(os.path.basename(path))[0].replace( '-benchmark', '' )

writer.writerow( map( shortname, tests ) )
writer.writerow( times )