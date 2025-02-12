# based on rename.tcl, translated to rename.py by wdecoster and copilot, 
# to be converted into an executable with pyinstaller --noconsole --onefile --noconfirm rename.py (in powershell)
# currently, that means for me: C:\users\wdecoster\appdata\local\packages\pythonsoftwarefoundation.python.3.9_qbz5n2kfra8p0\localcache\local-packages\python39\Scripts\pyinstaller.exe --onefile --console --noconfirm 'rename.py'
# the script works in the directory where it is located

import os
import glob

def main():
    print("This script helps you to change long filenames separated by underscores ( _ )")

    ext = input("Which files (extension) do you want to rename ?: ")
    filelist = glob.glob(f"*{ext}")
    lext = len(ext) + 1
    print(f"Found {len(filelist)} files with extension {ext}")
    print("Which parts of the original filename do you want ?")
  
    print(f"First file in the list ({filelist[0]}) has following parts:")
    parts = filelist[0].split('_')
    for i, part in enumerate(parts):
        print(f"part {i+1} : {part}")

    print("Type the numbers of the parts you want, separated by a space")
    outlist = list(map(int, input().split()))

    for file in filelist:
        lfile = len(file)
        i = lfile - lext
        out2 = file[:i]
        out = ""
        parts = out2.split('_')
        for ele in outlist:
            num = ele - 1
            breek = parts[num]
            out += breek + "_"
        out = out.rstrip('_')
        if os.path.exists(f"{out}.{ext}"):
            well = parts[1]
            out += f"_{well}"
        os.rename(file, f"{out}.{ext}")

if __name__ == "__main__":
    main()
