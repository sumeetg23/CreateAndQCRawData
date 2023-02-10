import os
import re


def validateSamplesheet(samplesheet):
    """
    Check for a few possible input errors and give an appropriate response
    """
    if not os.path.isfile(samplesheet):
        print("Couldn't be started. The following path or file %s does not exist" % (samplesheet,))
        return False
    return True


def validateOutput(outputfolder):
    """
    Check for a few possible input errors and give an appropriate response
    """
    if outputfolder.endswith('/'):
        outputfolder = outputfolder[:-1]
        if not os.access(dirname(outputfolder), os.W_OK):
            print("Couldn't be started. You don't have write permission for the following path %s" % (dirname(outputfolder),))
            print("It's also possible that parent folders/directories of the given path do not exist yet and need to be created")
            return None
    return outputfolder

