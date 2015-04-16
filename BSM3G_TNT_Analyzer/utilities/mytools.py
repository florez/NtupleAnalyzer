import ROOT as rt

def listRootFile(file,base=""):
    # the list
    path = []

    # get the base object
    base = base.strip("/")
    obj = file
    if len(base) != 0:
        obj = file.Get(base)
        if obj == None:
            return path
        path.append(base)
    
    # if base object is directory,
    # list recursively
    if obj.InheritsFrom("TDirectory"):
        for key in obj.GetListOfKeys():
            _base = (base + "/" + key.GetName()).lstrip("/")
            path.extend(listRootFile(file,_base))
    return path
