from fabric.api import env, local, run, sudo, cd, hide, prefix
from fabric.context_managers import shell_env, prefix
from fabric.operations import put, get
from fabric.contrib.files import append, exists
env.use_ssh_config = True

def all():
    prepare()
    numpy_release()

def prepare():
    prepare_apt()
    prepare_userspace()

def prepare_userspace():
    """
    This can be reverted by executing 'remove_userspace'.
    """
    copy()
    setup_wine()
    setup_paver()
    gitrepos()
    numpy_cpucaps()

def prepare_apt():
    #sudo("add-apt-repository -y ppa:ubuntu-wine/ppa") required for > precise install wine1.6
    sudo("apt-get -qq update")
    # This is needed to avoid the EULA dialog
    # (http://askubuntu.com/questions/16225/how-can-i-accept-the-agreement-in-a-terminal-like-for-ttf-mscorefonts-installer)
    sudo("echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections")
    sudo("apt-get -y install git wine patch python-virtualenv texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended make python-dev g++ libfreetype6-dev libpng-dev python-matplotlib python-setuptools python-jinja2 python-pygments python-docutils python-sphinx --no-install-recommends")
    # trusty cython, should probably be replaced with a regular upstream wheel
    sudo("dpkg -i /vagrant/cython_0.20.2-1_i386.deb")
    sudo("apt-get install -f")

def remove_userspace():
    """
    Deletes (!) the NumPy and Wine changes. Use with great care.
    """
    run("rm -rf repos")
    run("rm -rf .wine")

def copy():
    run("mkdir -p repos/numpy-vendor")
    run("cp -r /vagrant/* repos/numpy-vendor/")

def gitrepos():
    run("mkdir -p repos")
    with cd("repos"):
        run("git clone https://github.com/numpy/numpy")
        with cd("numpy"):
            run("git submodule init")
            run("git checkout -t origin/maintenance/1.10.x")
            run("git submodule update")

def setup_wine():
    with cd("repos/numpy-vendor"):
        run("./setup-wine.sh")

def setup_paver():
    with cd("repos/numpy-vendor"):
        run("tar xzf Paver-1.2.2.tar.gz")
        with cd("Paver-1.2.2"):
            run("python setup.py build")
            sudo("python setup.py install")

def numpy_cpucaps():
    with cd("repos/numpy/tools/win32build"):
            run('wine "C:\Python27\python" build-cpucaps.py')
            run(r"cp cpucaps.dll $HOME/.wine/drive_c/Program\ Files\ \(x86\)/NSIS/Plugins")

def numpy_release():
    with cd("repos/numpy"), shell_env(PYTHONPATH='/home/vagrant/repos/local/lib/python2.7/site-packages',
                                      NPY_SEPARATE_COMPILATION='0'):
        run("rm -rf ../local")
        run("paver sdist")
        run("python setup.py install --prefix ../local")
        run("paver pdf")
        run("paver bdist_superpack -p 3.4")
        run("paver bdist_superpack -p 3.3")
        run("paver bdist_superpack -p 2.7")
        run("paver write_release_and_log")
        run("paver bdist_wininst_simple -p 2.7")
        run("paver bdist_wininst_simple -p 3.3")
        run("paver bdist_wininst_simple -p 3.4")
    numpy_copy_release_files()

def numpy_copy_release_files():
    with cd("repos/numpy"):
        run("cp -r release/ /vagrant/")
        run("cp build_doc/pdf/*.pdf /vagrant/release/")
        run("mkdir /vagrant/release/pypi")
        run("cp dist/*.exe /vagrant/release/pypi/")
    with cd("/vagrant/release"):
        run("mv installers/* .")
        run("rm -r installers")
        run("mv NOTES.txt README.txt")

mac_tmp = "numpy_tmp" # NumPy and dependencies will be built in $HOME/"mac_tmp"
#FIXME: the Python path is hardwired here...
#TODO: just remove this hack:
mac_prefix = 'export PYTHONPATH="$HOME/%s/usr/lib/python2.7/site-packages/" PATH="$HOME/%s/usr/bin:$PATH"' % (mac_tmp, mac_tmp)

def mac_setup():
    run("mkdir %s" % mac_tmp)
    mac_setup_numpy()
    #mac_setup_paver()
    #mac_setup_virtualenv()
    #mac_setup_bdist_mpkg()
    mac_copy_pdf()

def mac_copy_pdf():
    with cd(mac_tmp + "/numpy"):
        run("mkdir -p build_doc/pdf")
        put("release/reference.pdf", "build_doc/pdf/")
        put("release/userguide.pdf", "build_doc/pdf/")

def mac_remove_userspace():
    run("rm -rf %s" % mac_tmp)

def mac_setup_numpy():
    with cd(mac_tmp):
        run("git clone https://github.com/numpy/numpy")
        with cd("numpy"):
            run("git submodule init")
            run("git checkout -t origin/maintenance/1.9.x")
            run("git submodule update")

def mac_setup_bdist_mpkg():
    with cd(mac_tmp):
        run("git clone https://github.com/rgommers/bdist_mpkg.git")
        with cd("bdist_mpkg"):
            with prefix(mac_prefix):
                run("python setup.py install --prefix=../usr")

def mac_setup_paver():
    with cd(mac_tmp):
        with prefix(mac_prefix):
            put("Paver-1.2.2.tar.gz", ".")
            run("tar xzf Paver-1.2.2.tar.gz")
            with cd("Paver-1.2.2"):
                run("python setup.py install --prefix=../usr")

def mac_setup_virtualenv():
    with cd(mac_tmp):
        with prefix(mac_prefix):
            put("virtualenv-1.8.4.tar.gz", ".")
            run("tar xzf virtualenv-1.8.4.tar.gz")
            with cd("virtualenv-1.8.4"):
                run("python setup.py install --prefix=../usr")

def mac_numpy_release():
    with cd(mac_tmp + "/numpy"):
        run("paver sdist")
        run("paver dmg -p 2.6")
        run("paver dmg -p 2.7")
        get("release/installers/*.dmg", "release/")

# ------------------------------------------------
# Vagrant related configuration

def vagrant():
    vc = _get_vagrant_config()
    # change from the default user to 'vagrant'
    env.user = vc['User']
    # connect to the port-forwarded ssh
    env.hosts = ['%s:%s' % (vc['HostName'], vc['Port'])]
    # use vagrant ssh key
    env.key_filename = vc['IdentityFile'].strip('"')
    # Forward the agent if specified:
    env.forward_agent = vc.get('ForwardAgent', 'no') == 'yes'

def _get_vagrant_config():
    """
    Parses vagrant configuration and returns it as dict of ssh parameters
    and their values
    """
    result = local('vagrant ssh-config', capture=True)
    conf = {}
    for line in iter(result.splitlines()):
        parts = line.split()
        conf[parts[0]] = ' '.join(parts[1:])
    return conf

# ---------------------------------------
# Just a simple testing command:

def uname():
    run('uname -a')
