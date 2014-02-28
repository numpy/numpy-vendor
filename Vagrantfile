# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant::Config.run do |config|
  #config.vm.box = "precise64"
  #config.vm.box_url = "http://files.vagrantup.com/precise64.box"
  config.vm.box = "precise32"
  config.vm.box_url = "https://cloud-images.ubuntu.com/vagrant/precise/current/precise-server-cloudimg-i386-vagrant-disk1.box"
# issues with shared folder
#config.vm.box = "saucy32"
#config.vm.box_url = "https://cloud-images.ubuntu.com/vagrant/saucy/current/saucy-server-cloudimg-i386-vagrant-disk1.box"
# needs newer wine (bug in 1.4.1), e.g. from ppa, cpuid does not build
#config.vm.box = "raring32"
#config.vm.box_url = "https://cloud-images.ubuntu.com/vagrant/raring/current/raring-server-cloudimg-i386-vagrant-disk1.box"
#config.vm.box = "trusty32"
#config.vm.box_url = "https://cloud-images.ubuntu.com/vagrant/trusty/current/trusty-server-cloudimg-i386-vagrant-disk1.box"
  config.ssh.forward_agent = true
end
