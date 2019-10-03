from invoke import task 
import os

class colors:
	BLUE = '\033[1;34;40m'
	GREEN =  '\033[32m'
	YELLOW =  '\033[33m'
	WHITE = '\033[37m'
	RESET = '\033[0;0m'

def print_color(color, msg):
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

	print_color(colors.GREEN, ">>> Cleaning! <<<")
	commands = [
		'cd ' + project_pwd,
		'rm -Rf build',
		'rm -Rf artifacts',
	]
	ctx.run(' && '.join(commands))
	print_color(colors.GREEN, ">>> build and artifact folders deleted. <<<")


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

	build_commands = [
		'cd ' + project_pwd,
		'mkdir build',
		'cd build',
		'mkdir makefiles',
		'cd makefiles',
		'cmake ../..',
		'cmake --build .',
	]
	if build_folder_exists:
		build_commands.remove('mkdir build')
		build_commands.remove('mkdir makefiles')

	print_color(colors.BLUE, ">>> Building! <<<" )
	ctx.run(' && '.join(build_commands))

	install_commands = [
		'cd ' + project_pwd,
		'cd build',
		'cd makefiles',
		'make install',
	]
	print_color(colors.BLUE, ">>> Installing! <<<" )
	ctx.run(' && '.join(install_commands))

	print_color(colors.BLUE, ">>> Builded and installed. <<<")


@task(
	help = {
		'verbose': "Print commands before execute them."
	}
)
def run_case_ex(ctx, verbose=False):
	"""
	Run an example case.
	"""
	project_name, project_pwd = get_project_name_and_folder()
	orig_folder = project_pwd + '/src/cpp/inCav.txt '
	dest_folder = project_pwd + '/artifacts/'

	dest_folder_exists = os.path.isdir(dest_folder) and os.path.exists(dest_folder)
	if not dest_folder_exists:
		build(ctx)
	
	commands = [
		'cp ' + orig_folder + dest_folder,
		'cd ' + dest_folder,
		'./' + project_name + 'Old inCav.txt',
	]

	if verbose:
		print(commands)

	print_color(colors.YELLOW, '>>> Running! <<<')
	ctx.run(' && '.join(commands))
	print_color(colors.YELLOW, '>>> Finished! <<<')
