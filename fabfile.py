from fabric.api import env, local, run, sudo, cd, hide
from fabric.contrib.files import append, exists

def prepare():
    sudo("apt-get -qq update")
    sudo("apt-get -qq install git")

def gitrepos():
    run("mkdir -p repos")
    with cd("repos"):
        run("git clone --depth=1 https://github.com/certik/numpy-vendor")

def setup_wine():
    with cd("repos/numpy-vendor"):
        run("./setup-wine.sh")

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
