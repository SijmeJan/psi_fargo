def read_file(file_name):
    f = open(file_name, "r")
    lines = f.readlines()
    f.close()

    return lines

def write_file(file_name, lines):
    f = open(file_name, "w")
    f.writelines(lines)
    f.close()

def add_line(lines, str_line):
    '''Add line to end of list'''
    n = len(lines)
    lines[n:n] = [str_line]
