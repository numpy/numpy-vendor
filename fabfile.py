from fabric.api import env, local, run, sudo, cd, hide
from fabric.contrib.files import append, exists

def all():
    prepare()
    numpy_release()

def prepare():
    prepare_apt()
    copy()
    setup_wine()
    setup_paver()
    gitrepos()
    numpy_cpucaps()

def prepare_apt():
    sudo("apt-get -qq update")
    # This is needed to avoid the EULA dialog
    # (http://askubuntu.com/questions/16225/how-can-i-accept-the-agreement-in-a-terminal-like-for-ttf-mscorefonts-installer)
    sudo("echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections")
    sudo("apt-get -qq install git wine python-virtualenv texlive-latex-recommended make")

def copy():
    run("mkdir -p repos/numpy-vendor")
    run("cp /vagrant/* repos/numpy-vendor/")

def gitrepos():
    run("mkdir -p repos")
    with cd("repos"):
        run("git clone https://github.com/numpy/numpy")
        with cd("numpy"):
            run("git checkout -t origin/maintenance/1.7.x")

def setup_wine():
    with cd("repos/numpy-vendor"):
        run("./setup-wine.sh")

def setup_paver():
    with cd("repos/numpy-vendor"):
        run("tar xzf Paver-1.0.5.tar.gz")
        with cd("Paver-1.0.5"):
            run("sudo python setup.py install")

def numpy_cpucaps():
    with cd("repos/numpy/tools/win32build/cpucaps"):
            run('wine "C:\Python27\python" "C:\Python27\Scripts\scons.py"')
            run(r"cp cpucaps.dll $HOME/.wine/drive_c/Program\ Files\ \(x86\)/NSIS/Plugins")

def numpy_release():
    with cd("repos/numpy"):
        run("time paver bdist_superpack -p 3.2")

# ------------------------------------------------
# Vagrant related configuration

def vagrant():
    vc = _get_vagrant_config()
    # change from the default user to 'vagrant'
    env.user = vc['User']
    # connect to the port-forwarded ssh
    env.hosts = ['%s:%s' % (vc['HostName'], vc['Port'])]
    # use vagrant ssh key
    env.key_filename = vc['IdentityFile']
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
