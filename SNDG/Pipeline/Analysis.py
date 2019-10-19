import logging
import traceback
from termcolor import colored

_log = loggin.getLogger('py.Analysis')

def execute(cmd_unformated,wd="./",retcodes=[0], **kargs):
    cmd = "cd " + wd + ";" + cmd_unformated.format(**kargs)
    _log.debug( cmd)
    try:
        sp.check_output(cmd, shell=True, stderr=sp.STDOUT)
        _log.debug(cmd + " -> OK")
    except CalledProcessError as ex:
        _log.warning(ex.output)
        if ex.returncode not in retcodes:
            raise

class Analysis():

    _default_input_folder = "input/"

    def __init__(self,name):
        self.name = name
        self.description = ""
        self.docker_image = None
        self.command = ""
        self.inputs = {}
        self.outputs = {}
        self.params = {}

    def execute(self,env):
        cb = CommandBuilder()
        cb.baseDir(env.base)

        for input_name,input_datasource in self.inputs.items():
            execute("ln {input} {output}",
                    output= Analysis._default_input_folder,
                    input=input_datasource.file_path())
            cb.bind_var(input_name,Analysis._default_input_folder + input_datasource.file_name())

        cmd = cb.build()

        execute(cmd)

        for input_name,input_datasource in self.inputs.items():
            execute("ln {input} {output}",
                    output= Analysis._default_input_folder,
                    input=input_datasource.file_path())
            cb.bind_var(input_name,Analysis._default_input_folder + input_datasource.file_name())





