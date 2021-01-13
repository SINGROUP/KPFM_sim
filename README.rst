================
KPFM simulation tools
================
Description
-----------

A set of Python modules for running DFT level KPFM simulations using CP2k and storing the results to a relational database (SQLite). Works for AFM simulations as well since they are just KPFM simulations without bias modulation. Contains three main components to help you with the simulations:
1. Planning CP2k calculation tasks and storing them to a database for later execution.
2. Running previously planned tasks automatically exploiting the results of the previous task step to reduce computation time.
3. Storing the essential results (total energy, atomic geometry...) in a structured form to a relational database and extracting those results easily within Python scripts.
These tools expect you to set up your atomic configuration using Atomic Simulation Environment (ASE). If you do not like ASE, you can use whatever modeling tools you like and read the atomic configuration with ASE to an Atoms object and use it.

Requirements
------------

- Python 3.x  (originally Python 2.7)
- Cython (`http://cython.org/ <http://cython.org/>`_)
- NumPy
- SciPy
- Atomic Simulation Environment (ASE), `https://wiki.fysik.dtu.dk/ase/ <https://wiki.fysik.dtu.dk/ase/>`_
- `CP2k_mtools <https://github.com/SINGROUP/CP2k_mtools> `_ (python base)
- `DFT_gridIO <https://github.com/SINGROUP/DFT_gridIO>`_ (python base)

For calculating electrostatic potential:
`KPFM_FEM_tools <https://github.com/SINGROUP/KPFM_FEM>`_ (python base)


WIKI
____

Newest description, installation guide and examples for running in the under development `wiki pages <https://github.com/SINGROUP/KPFM_sim/wiki>`_

Installation
------------

Put this directory containing the Python and Cython (.pyx) modules to your ``PYTHONPATH`` environment variable. The Cython modules should be compiled automatically at runtime if the line as long as you have Cython installed. If not, you can try running the setup scripts as
``python setup_x.py build_ext --inplace``
Create a new environment variable called ``KPFM_GLOBAL_SCRIPTS`` and store the path to the ``scripts`` folder to it. These are scripts that have no system dependent parameters so you can call them from a centralized location, and adding that location to an environment variable makes your life easier.

Usage
-----

This is an introduction to using these tools to run KPFM simulations. For more information about multiscale KPFM simulations in general, see the ``docs`` folder. If you have access to the shared archive folder of the SIN group, you can find a complete usage example in ``/path_to_archive/jritala/KPFM_sim``. You should use the scripts in ``control_script_templates`` folder as a starting point for your own simulations since they are part of the user interface to these tools. See also ``slurm_script_templates`` folder for examples of Slurm scripts to actually run the simulations on a supercomputer. The user interface to these tools is far from refined so I do not mind if you choose to improve it or write your own. The underlying functionality works very well, in my opinion.

Setting up initial atomic configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the ``kpfm_init_cu_tip_on_nacl.py`` script as a template for creating your initial atomic configuration. This configuration should contain a tip model at the maximum distance from the sample. It does not matter how you create the model, as long as you have it stored in an ASE Atoms object in the end. To optimize the initial geometry, call

::

	cp2k_initializer = CP2k_init(project_name, atoms)
	cp2k_calc = cp2k_initializer.init_desc_tip()
	cp2k_calc.run()

where ``project_name`` is whatever name you choose and ``atoms`` is an Atoms object containing your tip-sample model. See the documentation of `CP2k_tools <>`_ for more information. After that, you should read the optimized geometry and CP2k output file, label the "roles" of the atoms in the model and store all this information to a results database defined by ``Result_db`` class in the ``kpfm_sim_result_db`` module. Now you should have an SQLite database file containing information on atoms and unit cell of your model and the initial geometry in the directory where you ran the initialization script. ``Result_db`` class represents the database in Python and it contains many methods for extracting data from it without any knowlegde of relational databases. If you want to see the actual structure of the database file and its contents, use ``sqlite3`` command line tool (requires knowledge of SQL language) or ``Sqliteman`` which has a GUI.

Probing the sample at different positions (AFM simulation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now that you have initialized the result database and stored the initial atomic configuration to it, you can easily run a geometry optimization at multiple positions of the tip above the sample to get the total energy, atomic geometry and forces on atoms at those positions. These tools allow you to plan multiple tasks and store them to a database, and when you submit a batch job, one of these tasks is fetched for execution. Each task typically consists of multiple steps that each correspond to a single CP2k calculation. The task type to sample different positions above the sample is called "descend tip" and it is defined by the ``Descend_tip_task`` class in ``kpfm_sim_tasks`` module. The descend tip task samples a range of tip-sample distances starting from a position far away and moving the tip closer by constant amount at each step.

To create a descend tip task, use the ``plan_descend_tip_task.py`` script template in ``control_script_templates`` folder and change the parameters to match your simulation. Make sure ``V = 0`` to run without bias voltage. Then call

``python plan_descend_tip_task.py <task_db_file> <result_db_file> <global_res_db_file>``

where ``<task_db_file>``, ``<result_db_file>`` and ``<global_res_db_file>`` are relative paths to task, result and global result database files respectively.

**Important implementation note:**
In principle, the SQLite database handles concurrent writes to it correctly. However, in an environment with a parallel file system, as in many supercomputers like the CSC clusters, concurrent writes to the same database file from processes running on different nodes may result in corruption of the database. The workaround I have used is to create separate task and result database files for each concurrently running job and copy the results from the separate result databases to a global result database at the end. I suggest you create as many subfolders as there are concurrent jobs that you want to run, and call them ``worker_n``, for example, where ``n`` labels the different jobs.

If you want to execute multiple tasks at the same time in an environment with a parallel file system, the ``task_db_file`` and ``result_db_file`` should be separate for each concurrently running job as described in the implementation note above. If you follow the suggested scheme, you should call ``plan_descend_tip_task.py`` with arguments

``python plan_descend_tip_task.py worker_1/tasks.db worker_1/results.db your_simulation_results.db``

where ``worker_1`` is a subfolder you created and ``your_simulation_results.db`` is the database file containing the initial atomic configuration. ``tasks.db`` and ``results.db`` files are created automatically if they do not exist and the task you planned is saved to the ``tasks.db`` database file.

To execute a task you have planned, run the ``run_task.py`` script found in the ``scripts`` folder as

``python run_task.py <task_db_file> <project_path> <slurm_id> [type_constraint] [status_constraint]``

where ``<task_db_file>`` is a relative path from ``<project_path>`` to the task database file and ``<project_path>`` is the absolute path to the root directory of the simulations. ``project_path`` is needed in cases where the CP2k is run on a local file system of a node but the database files are on the shared file system. Since the tasks are typically executed in Slurm batch jobs, the ``<slurm_id>`` should be set to the ID of the slurm job executing the task. ``[type_constraint]`` and`` [status_constraint]`` are optional and can be used to restrict the type of the task to be run if there are multiple different kinds of tasks waiting and you want to run a specific one. See the ``worker_task_batch.sh`` script in ``slurm_script_templates`` for an example of a Slurm script (written for CSC Taito cluster). In particular, you should have the line
``trap "python $KPFM_GLOBAL_SCRIPTS/call_error_handler.py $SLURM_JOB_ID $ORIG_DIR $TASK_DB_FILE; exit" ERR TERM``
in the Slurm script if you want to have the error handler working. It is not necessary, but makes restarting possible in the case of an error or exceeded time limit. Otherwise you have to modify the task database by hand. You may have to do that anyway, if the cause of termination is something else than time limit. In that case, open the task database file using sqlite3 or Sqliteman and change the task state to waiting.

Probing the sample with different bias voltages (KPFM simulation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The way how the bias voltage between the probe and the sample holder is applied to the KPFM simulation depends on the type of the system you are studying. In particular, there are two entirely different cases:

1. Thick dielectric sample (thick meaning that you cannot model the whole sample within DFT)
2. Thin dielectric sample on metallic substrate (thin meaning that you can model the whole sample as well as some layers of the metal substrate)

In the case of a thick dielectric sample, you should calculate the electrostatic potential generated by the macroscopic part of the probe-sample model using `KPFM_FEM_tools <>`_. See the documentation of that package for instructions. When you have calculated the potential at sufficient range of tip-sample distances and have them stored into a FEM results database file, you should copy the data into the KPFM results database that was created during initialization of the atomic configuration. You can do that using ``copy_pot_to_result_db.py`` script. When you execute a task with a non-zero bias voltage, the electrostatic potential is read from the database and written into a cube file by a function in the ``axisym_pot_to_cube`` module. That cube file is read by CP2k and added as an external potential to the DFT calculation.

If you have a thin sample, however, the electrostatic potential between the tip and the metallic substrate is entirely defined by the atomic model. The correct potential/field between the tip and the substrate is generated if a suitable amount of charge is transferred between them. This happens if one is able to shift the Fermi levels of the tip and the substrate with respect to each other. One way to do this is to apply step-like external potential to the DFT calculation so that the tip is at a different potential than the substrate. You can use ``piecewise_linear_potential`` module to create the step-like potential. There is no option to use this method automatically within the simulation tools environment yet. Find the TODO comment in ``kpfm_sim_tasks`` if you want to implement it.

Independent of the way of applying the bias voltage, you can either go through the zero bias scan points and vary the bias at each of those points or fix the bias voltage and descend the tip with that bias. You can use the ``plan_tune_bias_task.py`` or ``plan_tune_bias_tasks_srange.py`` as a template for a script for planning tasks which have varying bias voltage. Descending the tip using a fixed bias voltage works by planning tasks using ``plan_descend_tip_task.py`` script with a non-zero ``V``. Descending with a fixed bias seems to work better because varying the bias changes the atomic geometry globally and thus the previous step with a different bias voltage is not a good guess for the initial geometry of a geometry optimization.

Calculating atomic forces
^^^^^^^^^^^^^^^^^^^^^^^^^
The forces on atoms must be calculated during a separate run because the forces on fixed atoms are zero during geometry optimization. Use ``calc_atomic_forces.py`` script to do it.

Combining the results into one database file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you executed multiple tasks in parallel and have multiple separate database files, you can combine them into one database using the ``copy_scan_points.py`` script.

Analysing the results in the database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Result_db class defined in kpfm_sim_result_db module contains many methods for extracting data from the SQLite result database without any knowledge of relational databases. You can also use the ready-made ``extract_*`` scripts in the scripts folder or use them as an example.

Author
------
Ondrej Krejci (2021)
`ondrej.krejci@aalto.fi <mailto:ondrej.krejci@aalto.fi>`

Juha Ritala (2016)
`jritala@gmail.com <mailto:jritala@gmail.com>`_

