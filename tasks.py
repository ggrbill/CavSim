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


def check_option(option, option_list):
    """
    :return: Return True if an option is valid in a given option list
    """
    if option not in option_list:
        print_color(colors.YELLOW, f' The option \'{option}\' is not valid. It should be one of {option_list}.')
        return False

    return True


@task()
def clean(ctx):
	"""
	Delete 'build' and 'artifacts' folders.
	"""
	current_dir = os.path.abspath(os.curdir)
	if 'build' in current_dir or 'artifacts' in current_dir:
		print_color(colors.YELLOW, 'Impossible to delete build or artifacts folder. Your current directory is one of them.')
		return

	project_name, project_pwd = get_project_name_and_folder()

	print_color(colors.GREEN, ">>> Cleaning! <<<")
	commands = [
		'cd ' + project_pwd,
		'rm -Rf build/_ninja',
		'rm -Rf build/_makefile',
		'rm -Rf artifacts',
	]
	ctx.run(' && '.join(commands))
	print_color(colors.GREEN, ">>> build and artifact folders deleted. <<<")


@task(
	help = {
		'cclean': "Call 'clean' task (Delete 'build' and 'artifacts' folders) before build again.",
		'sys': "The build system can be 'makefile' or 'ninja'(default).",
	}
)
def build(ctx, cclean=False, sys='ninja'):
	"""
	Build C++ code and install the artifacts.
	"""
	if not check_option(sys, ['makefile', 'ninja']):
		return

	sys_build = {
		'makefile' : {
			'Generate' : '-G"Unix Makefiles"',
			'Install' : 'make install',
		},
		'ninja' : {
			'Generate' : '-G"Ninja"',
			'Install' : 'ninja install',
		},
	}

	project_name, project_pwd = get_project_name_and_folder()

	if cclean:
		clean(ctx)
	
	build_folder = project_pwd + '/build/'
	build_folder_exists = os.path.isdir(build_folder) and os.path.exists(build_folder)
	
	sys_build_folder = project_pwd + '/build/_' + sys
	sys_build_folder_exists = os.path.isdir(sys_build_folder) and os.path.exists(sys_build_folder)

	build_commands = [
		'cd ' + project_pwd,
		'mkdir build',
		'cd build',
		'mkdir _' + sys,
		'cd _' + sys,
		'cmake ' + sys_build[sys]['Generate'] + ' ../..',
		'cmake --build .',
	]
	if build_folder_exists:
		build_commands.remove('mkdir build')
	if sys_build_folder_exists:
		build_commands.remove('mkdir _' + sys)
	print_color(colors.BLUE, ">>> Building! <<<" )
	ctx.run(' && '.join(build_commands))
	
	install_commands = [
		'cd ' + project_pwd,
		'cd build',
		'cd _' + sys,
		sys_build[sys]['Install']
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
