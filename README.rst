How To Use
==========

Do::

    vagrant up
    fab vagrant prepare

Log in and manually inspect the numpy repository (do some last minute fixes,
etc.)::

    vagrant ssh
    cd repos/numpy
    # Do any changes that are not in official repositories

Do the release (build general and windows binaries)::

    fab vagrant numpy_release

And the directories ``release``, ``dist`` and ``build_doc`` will be copied to
the current directory from the VM. If you need anything else, just login using
``vagrant ssh`` and copy it to ``/vagrant`` inside the VM.
