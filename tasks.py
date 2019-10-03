from invoke import task 
import os

class colors:
	BLUE = '\033[1;34;40m'
	GREEN =  '\033[32m'
	YELLOW =  '\033[33m'
	WHITE = '\033[37m'
	RESET = '\033[0;0m'

def print_(color, msg):
	print(color + msg + colors.RESET)

def get_project_name_and_folder():
    """
    :return: The Name and Root Directory of the current project from the current working dir.
    """
    filename = 'environment.devenv.yml'
    directory = os.path.abspath(os.curdir)
    while True:
        if os.path.exists(os.path.join(directory, filename)):
            return os.path.basename(directory), directory.replace('\\', '/')

        directory, current_dir_name = os.path.split(directory)
        if len(current_dir_name) == 0:
            raise RuntimeError(
                'Can not obtain a project from the current working directory')
        assert len(directory) != 0


@task()
def clean(ctx):
	"""
	Delete 'build' and 'artifacts' folders.
	"""
	project_name, project_pwd = get_project_name_and_folder()
	ctx.run('cd ' + project_pwd)

	print_(colors.GREEN, ">>> Cleaning! <<<")
	commands = [
		'cd ' + project_pwd,
		'rm -Rf build',
		'rm -Rf artifacts',
	]
	ctx.run(' && '.join(commands))
	print_(colors.GREEN, ">>> build and artifact folders deleted. <<<")


@task(
	help = {
		'cclean': "Call 'clean' task (Delete 'build' and 'artifacts' folders) before build again."
	}
)
def build(ctx, cclean=False):
	"""
	Build C++ code and install the artifacts.
	"""
	project_name, project_pwd = get_project_name_and_folder()

	if cclean:
		clean(ctx)
	
	build_folder = project_pwd + '/build' 
	build_folder_exists = os.path.isdir(build_folder) and os.path.exists(build_folder)

	commands = [
		'cd ' + project_pwd,
		'mkdir build',
		'cd build',
		'mkdir makefiles',
		'cd makefiles',
		'cmake ../..',
		'cmake --build .',
		'make install',
	]
	if build_folder_exists:
		commands.remove('mkdir build')
		commands.remove('mkdir makefiles')

	print_(colors.BLUE, ">>> Building And Installing! <<<" )
	ctx.run(' && '.join(commands))
	print_(colors.BLUE, ">>> Builded and installed. <<<")


@task()
def run_case_ex(ctx):
	"""
	Run an example case.
	"""
	project_name, project_pwd = get_project_name_and_folder()
	orig_folder = project_pwd + '/src/cpp/inCav.txt '
	dest_folder = project_pwd + '/artifacts/'

	commands = [
		'cp ' + orig_folder + dest_folder,
		'cd ' + dest_folder,
		'./' + project_name + 'Old inCav.txt',
	]
	print_(colors.YELLOW, '>>> Running! <<<')
	print(commands)
	ctx.run(' && '.join(commands))
	print_(colors.YELLOW, '>>> Finished! <<<')
