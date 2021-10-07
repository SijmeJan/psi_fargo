def read_file(file_name):
    '''Return all lines in text file'''
    f = open(file_name, "r")
    lines = f.readlines()
    f.close()

    return lines

def write_file(file_name, lines):
    '''Write a list of strings to text file'''
    f = open(file_name, "w")
    f.writelines(lines)
    f.close()

def add_line(lines, str_line):
    '''Add line to end of list'''
    n = len(lines)
    lines[n:n] = [str_line]
