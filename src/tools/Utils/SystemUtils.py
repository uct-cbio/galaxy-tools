#! /usr/bin/python
import sys, os

def find_program(program):
    found_path = False
    program_path = ""
    for path in os.getenv("PATH").split(":"):
        check_path =  path + "/" + program
        if (os.path.isfile(check_path)):
            program_path = check_path
            found_path = True
            break
    if found_path:
        return program_path
    else:
        print "Can't find the " + program + " utility on your system"
        print "Please ensure it is installed and in your path"
        sys.exit()

if __name__ == "__main__":
    program = find_program(sys.argv[1])
    print program


