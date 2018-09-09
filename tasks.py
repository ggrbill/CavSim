from invoke import task 
import os

def get_project_name_and_folder():
    """
    :return: The Name and root directory of the current project from the current working dir.
    """
    filename = 'environment.yml'
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
	Delete 'build' folder.
	"""
	project_name, project_pwd = get_project_name_and_folder()
	ctx.run('cd ' + project_pwd)

	print("Cleaning!")
	commands = [
		'cd ' + project_pwd,
		'rm -Rf build',
	]
	ctx.run(' && '.join(commands))


@task(
	help = {
		'cclean': "Call Clean task (Delete 'build' folder) before build again."
	}
)
def build(ctx, cclean=False):
	"""
	Build C++ code.
	"""
	if cclean:
		clean(ctx)
	
	project_name, project_pwd = get_project_name_and_folder()
	
	print("Building!")
	commands = [
		'cd ' + project_pwd,
		'mkdir build',
		'cd build',
		'mkdir makefiles',
		'cd makefiles',
		'cmake ../..',
		'cmake --build .',
	]
	ctx.run(' && '.join(commands))


@task()
def run_case_ex(ctx):
	"""
	Run an example case.
	"""
	project_name, project_pwd = get_project_name_and_folder()
	orig_folder = project_pwd + '/src/cpp/inCav.txt '
	dest_folder = project_pwd + '/build/makefiles/bin '

	commands = [
		'cp ' + orig_folder + dest_folder,
		'cd ' + dest_folder,
		'./' + project_name + ' inCav.txt',
	]
	print(commands)
	ctx.run(' && '.join(commands))
