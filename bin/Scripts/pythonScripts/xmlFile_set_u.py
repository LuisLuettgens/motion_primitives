import sys
import re
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove, getcwd

def replace(file_path, uni_line_identifier,pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'r+') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                if uni_line_identifier in  line:
                    new_file.write(line.replace(pattern, subst))
                else:
                    new_file.write(line)
                
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)


if __name__ == "__main__":
    
    filename         = sys.argv[1]
    velocity         = sys.argv[2]
    unique_line_str  = sys.argv[3]
    bound            = sys.argv[4]

    bound_val = ""    
    
    with open(filename, "r+") as inf:
        for line in inf:    
            if unique_line_str in  line:
                line_local = line
                splitted=re.split(' ',line)
                for splinter in splitted:
                    if bound+'=' in splinter:
                        bound_split=re.split('=',splinter)
                        bound_val = bound_split[1].replace("\"", "")
                        
                            
        
    replace(filename,unique_line_str, bound +'='+ '"'+bound_val+ '"', bound+'='+ '"' + velocity + '"')
    print bound_val
