
=======
Running
=======

The reproductive simulation library, as it's name suggests, is a library and as such the concept of running the library doesn't really exist.  However you can make use of the applications that rely on this library.  A small repository of applications is available from `here <https://github.com/VirtualPregnancy/reprosim-examples>`_.  

Making use of virtual environments
==================================

Possibly the best way to make use of the Reprosim library is through the python bindings made available through virtual environments.  This allows a correspondance of library configuration to virtual environment that can be easily moved between i.e. for each build configuration Debug, Release we create a corresponding virtual environment.

Follow these steps for creating a Python virtual environment from which the Reprosim library will be available.

Create a home for all virtual environments
------------------------------------------

The first task is to create a directory to hold the virtual environment installations::

  mkdir virtual_environments
  
This directory can be created anywhere.

Create a virtual environment
----------------------------

The second task is to create a Python virtual environment to install the Reprosim python modules into::

  cd virtual_environments # change directory to where the virtual environment should be created
  virtualenv --system-site-packages develop
  
The *--system-site-packages* flag allows the virtual environment to access all the packages that the system python has installed.  This is useful for big packages which may be required, for example; numpy or scipy.  The name of the virtual environment (in this case *develop*) is determined from the branch of the Reprosim library that is going to be available.

Activate virtual environment
----------------------------

The third task is to activate the Python environment.  This can be done by executing a shell script made available in the installation, for POSIX systems execute the command::

  source /path/to/env/bin/activate
  
for Windows the equivalent command is::

  \path\to\env\Scripts\activate
  
The activate script may alter the command prompt to indicate the active virtual environment.  This script will also make changes to your path variables.  To undo these changes execute the *deactivate* script::

  deactivate
  
Install Reprosim into virtual environment
-----------------------------------------

With an active Python virtual environment change directory into the reprosim build directory::

  cd /path/to/reprosim-build/
  
From this directory change into the *src/bindings/python* directory::

  cd src/bindings/python
  
in this directory a Python file named *setup.py* should exist.  To make the Reprosim library available via the active virtual environment execute the following command::

  python setup.py develop
  
This will create a link from the active virtual environment to the Reprosim library, thus making the Reprosim python library available from the the currently active Python environment.

Test Reprosim in virtual environment
------------------------------------

With the virtual environment active that Reprosim is linked to, run python to get a command prompt::
  
  python
  
and at the command prompt enter the following::

  >>> from reprosim.diagnostics import set_diagnostics_level
  
if all has gone correctly ... nothing should happen! Another command prompt should appear::

  >>>

If the above command was successful then the Python applications given above will also run successfully.

Finally
-------

This procedure of making the Reprosim library available through a python virtual environment can be repeated for different builds of the Reprosim library.  The virtual environment is lightweight and provides great encapsulation for development of the library.
