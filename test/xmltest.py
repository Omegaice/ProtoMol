import glob
import os
import subprocess
import comparator
import optparse
import logging

DEFAULT_EPSILON = 0.00001
DEFAULT_SCALINGFACTOR = 1.0

def parse_params(flname):
    fl = open(flname)
    params = dict()
    for ln in fl:

        if ln[:2] != '##':
            continue
        eq_pos = ln.find('=')
        if eq_pos == -1:
            continue
        param_name = ln[2:eq_pos].strip()
        param_str = ln[eq_pos + 1:].strip()
        param_val = eval(param_str)
        params[param_name] = param_val
    fl.close()
    return params

# Arguments
parser = optparse.OptionParser(description='ProtoMol Test Suite')
parser.add_option('--verbose', '-v', action='store_true', dest='verbose', default=False, help='Verbose output')
parser.add_option('--parallel', '-p', action='store_true', dest='parallel', default=False, help='MPI Testing')

(options, args) = parser.parse_args()

# Setup Statistics
stats_test = 0
stats_pass = 0
stats_fail = 0
stats_error = 0

# Empty Output Directory
outdata = glob.glob( os.path.join(os.getcwd(),'tests','output','*' ) )
map( os.remove, outdata )

# Find Executable
executable = []
if options.parallel:
    executable.append( "mpirun" )
    executable.append( "-np" )
    executable.append( "2" )
executable.append( os.path.join(os.getcwd(), 'ProtoMol') )

# Verbose
if options.verbose:
    logging.basicConfig(level=logging.DEBUG)

# Find Tests
tests = glob.glob(os.path.join(os.getcwd(), 'tests', '*.conf'))
tests.sort()

for test in tests:
	if test.find( "-benchmark" ) != -1:
		tests.remove( test )

stats_test = len(tests)

# Run Tests
testid = 0
failed_tests = []
passed_tests = []
disabled_tests = []

for test in tests:
    conf_param_overrides = parse_params(test)
    epsilon = conf_param_overrides.get('epsilon', DEFAULT_EPSILON)
    scaling_factor = conf_param_overrides.get('scaling_factor',
        DEFAULT_SCALINGFACTOR)

    testid += 1
    name = os.path.splitext(os.path.basename(test))[0]

    # Run ProtoMol
    cmd = executable[:]
    cmd.append( test )

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()

    if p.returncode > 0:
        stats_error += 1
        disabled_tests.append( { 'id': testid, 'name': name }  )
        continue

    # Find Outputs
    expects = []
    outputs = glob.glob('tests/output/' + name + '.*')
    outputtemp = []
    for output in outputs:
        outbase = os.path.basename(output)
        if os.path.exists('tests/expected/' + outbase):
            outputtemp.append(output)
            expects.append('tests/expected/' + outbase)

    outputs = outputtemp[:]

    # Remove Untestable Types
    def untestable(path):
        ftype = os.path.splitext(os.path.basename(path))[1]
        if ftype in ['.dcd', '.header', '.xtc']:
            return False
        else:
            return True

    expects = filter(untestable, expects)
    outputs = filter(untestable, outputs)

    # Compare Outputs
    output_size = len(outputs)
    output_pass = []
    output_fail = []

    for i in xrange(output_size):
        ignoreSign = False
        ftype = os.path.splitext(os.path.basename(outputs[i]))[1]

        if ftype == '.vec':
            ignoreSign = True

        if comparator.compare(expects[i], outputs[i], epsilon, scaling_factor, ignoreSign):
            string = outputs[i] + " matches"
            output_pass.append( string )
            print string
        else:
            string = outputs[i] + " differs"
            output_fail.append( string )
            print string

    print name + " " + str(output_size) + " " + str(len(output_pass))

    # Create XML
    if len(output_pass) == output_size:
        stats_pass += 1
        passed_tests.append( { 'id': testid, 'name': name } )
    else:
        stats_fail += 1

        def combine(a,b):
            return a + "\n" + b

        message = reduce(combine, output_fail)
        failed_tests.append( { 'id': testid, 'name': name, 'message': message } )

print ""

# Print Information
if len( failed_tests ) > 0:
    print "Failed: "

    for test in failed_tests:
        print str(test['id']) + ":" + test['name']

if len(disabled_tests) > 0:
    print "Disabled:"

    for test in disabled_tests:
        print str(test['id']) + ":" + test['name']

# Write XML File
fXML = open( "results.xml", "w" )
fXML.write("<?xml version=\"1.0\" encoding='ISO-8859-1' standalone='yes'?>\n")
fXML.write("<TestRun>\n")

fXML.write("\t<FailedTests>\n")
for test in failed_tests:
    fXML.write("\t\t<FailedTest id=\"" + str(test['id']) + "\">\n")
    fXML.write("\t\t\t<Name>" + test['name'] + "</Name>\n")
    fXML.write("\t\t\t<FailureType>Assertion</FailureType>\n")
    fXML.write("\t\t\t<Message>" + test['message'] + "</Message>\n")
    fXML.write("\t\t</FailedTest>\n")
fXML.write("\t</FailedTests>\n")

fXML.write("\t<SuccessfulTests>\n")
for test in passed_tests:
    fXML.write("\t\t<Test id=\"" + str(test['id']) + "\">\n")
    fXML.write("\t\t\t<Name>" + test['name'] + "</Name>\n")
    fXML.write("\t\t</Test>\n")
fXML.write("\t</SuccessfulTests>\n")

fXML.write("\t<Statistics>\n")
fXML.write("\t\t<Tests>" + str(stats_test) + "</Tests>\n")
fXML.write("\t\t<FailuresTotal>" + str(stats_fail) + "</FailuresTotal>\n")
fXML.write("\t\t<Errors>" + str(stats_error) + "</Errors>\n")
fXML.write("\t\t<Failures>" + str(stats_fail) + "</Failures>\n")
fXML.write("\t</Statistics>\n")

fXML.write("</TestRun>\n")

fXML.close()
