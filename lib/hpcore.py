#############################################################################
# Support functions for the HyPhy parser.
# 11.2020
#############################################################################

import sys, datetime, random, string
import timeit

#############################################################################

def readMeta(metafile):
    features, first = {}, True;
    for line in open(metafile):
        if line[0] == "#":
            continue;
        if first:
            first = False;
            continue;
        line = line.strip().split("\t");
        fid, ttype, chrome, start, end, strand, gid, num_cds = line;
        features[fid] = { 'gid' : gid, 'chrome' : chrome, 'start' : start, 'end' : end, 'strand' : strand, 'num-cds' : int(num_cds) };
    return features;

#############################################################################

def runTime(msg=False, writeout=False):
	if msg:
		if not msg.startswith("#"):
			msg = "# " + msg;
		PWS(msg, writeout);

	PWS("# PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])), writeout)
	PWS("# Script call:    " + " ".join(sys.argv), writeout)
	PWS("# Runtime:        " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"), writeout);
	PWS("# ----------------", writeout);

#############################################################################

def PWS(o_line, o_stream=False, std_stream=True):
# Function to print a string AND write it to the file.
	if std_stream:
		print(o_line);
	if o_stream:
		o_stream.write(o_line + "\n");

#############################################################################

def printWrite(o_name, v, o_line1, o_line2="", pad=0):
# Function to print a string AND write it to the file.
    if o_line2 == "":
        outline = o_line1;
    else:
        outline = o_line1 + " "*(pad-len(o_line1)) + o_line2;
    if v in [-1,1,2]:
        print(outline);
    if v in [-1,1]:
        f = open(o_name, "a");
        f.write(outline + "\n");
        f.close();

#############################################################################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a string to make it a given length
	spaces = sep * (totlen - len(string));
	return string + spaces;

#############################################################################

def getDate():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y");

#############################################################################

def getTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%H:%M:%S");

#############################################################################

def getDateTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

#############################################################################

def getRandStr(strlen=6):
# This function generates a random string to add onto the end of tmp files to avoid possible overwrites.
        return ''.join(random.choice(string.ascii_letters) for m in range(strlen));

#############################################################################

def report_step(step, step_start_time, step_status, prog_start_time, start=False, full_update=False):
# Uses psutil to gather memory and time info between steps and print them to the screen.

    globs = {};
    globs['psutil'] = False;

    dashes = 150
    if globs['psutil']:
        import psutil;
        dashes = 175;
    # Determine the number of dashes to frame the update table depending on the presence of psutil

    cur_time = timeit.default_timer();
    # The time at the start of the status update

    col_widths = [ 14, 10, 40, 40, 20, 16 ];
    if globs['psutil']:
        col_widths += [25, 20];
    # The column widths

    if start:
        headers = [ "# Date", "Time", "Current step", "Status", "Elapsed time (s)", "Step time (s)" ];
        if globs['psutil']:
            headers += ["Current mem usage (MB)", "Virtual mem usage (MB)"]
        # A list of the headers

        headers = "".join([ spacedOut(str(headers[i]), col_widths[i]) for i in range(len(headers)) ]);
        # Converting the list to a string based on the column widths

        printWrite("", 2, "# " + "-" * dashes);
        printWrite("", 2, headers);
        printWrite("", 2, "# " + "-" * dashes);
        # Print the dashes and the headers
    # The first call is just to print the headers

    ##########

    else:
        prog_elapsed = str(round(cur_time - prog_start_time, 5));
        # Get the total amount of time that the program has been running

        if not step_start_time:
        # If no step start time is given, then this is the first entry for this status
        # update, that will display "In progress..." or similar.

            out_line = [ "# " + getDate(), getTime(), step, step_status ];
            # The output for the initial status entry includes the date, time, step label, and progress message

            term_col_widths = col_widths[:4];
            # Only get the first 4 column widths for the initial status entry

            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];

            if full_update:
                out_line += "\n";
            # For some status updates, intermediate info will be printed, in which case we add a newline here

            sys.stdout.write("".join(out_line));
            sys.stdout.flush();
            # Convert the output list to a string, write, and flush stdout

        # The initial status entry to display "In progress..."

        #####

        else:
            step_elapsed = str(round(cur_time - step_start_time, 5));
            # Get the full step time here

            out_line = [ step_status, prog_elapsed, step_elapsed ];
            # Gather info for the full output line to print to screen

            if globs['psutil']:
                mem = round(sum([p.memory_info()[0] for p in globs['pids']]) / float(2 ** 20), 5);
                vmem = round(sum([p.memory_info()[1] for p in globs['pids']]) / float(2 ** 20), 5);
                out_line += [str(mem), str(vmem)];
            # If psutil is present, get current memory info

            term_col_widths = col_widths[3:];
            # Get the column widths for the print to screen output

            file_line = [ "# " + getDate(), getTime(), step ] + out_line;
            file_col_widths = col_widths[:3] + [30] + col_widths[4:];
            # For output to the file, we write the whole line each time
            # Add the initial entry fields here
            # This will also be used for some status updates where the whole message needs to be printed
            # to the screen
            
            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];
            file_line = [ spacedOut(str(file_line[i]), col_widths[i]) for i in range(len(file_line)) ];
            # Compile both the truncated and the full status update

            if full_update:
                sys.stdout.write("".join(file_line) + "\n");
                sys.stdout.flush();
            else:         
                sys.stdout.write("\b" * 40);
                sys.stdout.write("".join(out_line) + "\n");
                sys.stdout.flush();
            # For full updates, print the full line to the screen
            # For others, delete the "In progress..." column and update the same status line
            
            #printWrite(globs['logfilename'], 3, "".join(file_line));
            # Write the full line to the file.

        # The final status entry

        #####

    return cur_time;

#############################################################################

        # if scaff == "ScmyWZ3_7747_HRSCAF_7900":
        #     if int(end) < 35225145:
        #         scaff = "ScmyWZ3_7747_HRSCAF_7900_R";
        #         chrome = "chrX_R";
        #     elif int(start) >= 35225145:
        #         scaff = "ScmyWZ3_7747_HRSCAF_7900_NR"
        #         chrome = "chrX_NR";
        #     else:
        #         chrome = "chrX";
        # elif scaff in scaff_to_chr:
        #     chrome =  scaff_to_chr[scaff];
        # else:
        #     chrome = "NA";
