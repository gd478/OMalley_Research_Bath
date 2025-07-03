from MDANSE import REGISTRY

#Run in MDANSE command shell, likely found in downloads/Windows_executable_nightly/MDANSE\MDANSE_command_shell#

#Don't change anything not commented#

################################################################
# Job parameters                                               #
################################################################

parameters = {}
parameters['axis_selection'] = '1' #You'll need to open mdanse gui, open your trajectory file as .nc, open 'analysis > dynamics > Angular Correlation Function'#
                                   #select new definition and the atoms you want. Call the definition 1 2 3 etc for best compatibility with this template script.#
                                   
parameters['frames'] = (0, 999, 1) #(starting frame, final frame, step size)#
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\1', (u'hdf',)) #Change to your file path
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc' #Change to your file path

################################################################
# Setup and run the analysis                                   #
################################################################

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True) 

#Second AC job starts below, I have it calculating the acf for another axis over the same frames 0-999. AC jobs separated by '#' below#

parameters = {}
parameters['axis_selection'] = '2'
parameters['frames'] = (0, 999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\2', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '3'
parameters['frames'] = (0, 999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\3', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '1'
parameters['frames'] = (1000, 1999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\4', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '2'
parameters['frames'] = (1000, 1999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\5', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '3'
parameters['frames'] = (1000, 1999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\6', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '1'
parameters['frames'] = (2000, 2999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\7', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '2'
parameters['frames'] = (2000, 2999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\8', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '3'
parameters['frames'] = (2000, 2999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\9', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '1'
parameters['frames'] = (3000, 3999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\10', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '2'
parameters['frames'] = (3000, 3999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\11', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)

#

parameters = {}
parameters['axis_selection'] = '3'
parameters['frames'] = (3000, 3999, 1)
parameters['output_files'] = (u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\rot\\12', (u'hdf',))
parameters['per_axis'] = False
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'C:\\Users\\gd478\\Documents\\phd_yr2\\MD\\For_ILL+PSI\\building\\loaded\\bea\\2%\\20ps_post_prod\\bea2%.nc'

ac = REGISTRY['job']['ac']()
ac.run(parameters,status=True)