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
