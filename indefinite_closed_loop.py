# 
# Indefinite closed loop
#
# Authors:      Peter Weir and Marie Suver
# Last updated: September 1 2009

import sys
import timeimport numpy
import comedi
import PControl.Panel_com as Panel_com
from simple_step import Simple_Step
import curses

# -------------------------------------------------------------------------
# --------- experimenter input
# -------------------------------------------------------------------------
userInput = sys.argv

if len(userInput) < 2:
	print 'user input should be [UV or gr] [pol, qwp, qtr] [none]'
	sys.exit()

wavelength = userInput[1]   # UV ('UV') or green ('gr') LED 
filterType = userInput[2]   # 'pol' or 'qwp' or 'qtr'
eyeTreatment = userInput[3] # Will do no eye treatments (i.e. waxing) in this experiment, but will indicate 'none' anyways
	
if (wavelength != 'UV' and wavelength != 'gr') or (filterType != 'pol' and filterType != 'qwp' and filterType != 'qtr') or (eyeTreatment != 'none'):
	print 'user input should be [UV or gr] [qwp or qtr or pol] [none]'
	sys.exit()

# -------------------------------------------------------------------------
# --------- curses input initialization
# -------------------------------------------------------------------------
stdscr = curses.initscr()
#curses.nocbreak()
stdscr.keypad(0)
curses.echo()

stdscr.nodelay(1)
gin = stdscr.getch()

# -------------------------------------------------------------------------
# --------- DAQ initialization
# -------------------------------------------------------------------------
chans = [0,1]
subdev = 0
chrange = 0


# Open and configure device
daq = comedi.comedi_open('/dev/comedi0')
if not daq:
    raise "Error openning Comedi device"

maxdata = [0,0]
cr = [0,0]
maxdata[0] = comedi.comedi_get_maxdata(daq,subdev, chans[0])
cr[0]      = comedi.comedi_get_range(daq,subdev,chans[0],chrange)
maxdata[1] = comedi.comedi_get_maxdata(daq,subdev, chans[1])
cr[1]      = comedi.comedi_get_range(daq,subdev,chans[1],chrange)

# -------------------------------------------------------------------------
# --------- polarizer initialization
# -------------------------------------------------------------------------
NUM_TEETH_PASSIVE_GEAR = 80.0 # number of teeth on gear surrounding polarizer
NUM_TEETH_DRIVE_GEAR   = 64.0   # number of teeth on gear attached to motor shaft
STEPS_PER_ROT_MOTOR    = 6400.0  # setable on driver
RAMP_ACCEL             = 15000

stepsPerRotation = STEPS_PER_ROT_MOTOR*NUM_TEETH_PASSIVE_GEAR/NUM_TEETH_DRIVE_GEAR

# Open device 
dev = Simple_Step()

# enable device
dev.enable()

# No zeroing of device! Position calibrated (manually and occasionally by experimenters).

# DIO pin connected to LED
ledPin = 2
dev.set_dio_lo(ledPin) # turn LED off


# -------------------------------------------------------------------------
# --------- panels initialization
# -------------------------------------------------------------------------
NUM_PIXELS = 18*8
degreesPerPixel = 360.0/NUM_PIXELS

userport = '/dev/ttyS0'
panels = Panel_com.Panel_com( userport )
patternID = 1
panels.SetPatternID(patternID)

gainx = 16
gainy = 0
offsety = 0
offsetx = 0
modex = 1   # closed loop
modey = 0
patternYPos = 0

# -------------------------------------------------------------------------
# --------- closed loop
# -------------------------------------------------------------------------

direction = ['negative', 'positive']

TOTAL_TIME_VISUAL   = 10    # seconds of closed loop stimulus
TOTAL_TIME_EXP      = 48000 # seconds of closed loop (experimental) stimulus (many many minutes)
PRE_LED_TIME        = 5     # seconds LED off before turning on
PRE_MOTION_TIME     = 5     # seconds LED on before starting closed loop
trial_time = 0

stdscr.addstr('visual closed loop\n')
stdscr.refresh()
dev.set_dio_lo(ledPin)      # turn LED off
panels.SetPositions(143, 0)  # load [single] stripe pattern
panels.SetMode(modex, modey)
panels.SetGainOffset(gainx,offsetx,gainy,offsety)
panels.Start()
time.sleep(TOTAL_TIME_VISUAL)
panels.Stop()
panels.AllOff()


gain = -2500

start_pos_deg = numpy.random.randint(360)
start_pos = start_pos_deg*stepsPerRotation/360

date_time = time.strftime("%Y%m%d_%H%M%S")
fd = open('./' + wavelength + '_' + filterType + '_' + eyeTreatment + '_' + date_time + '.txt',mode='w')
fd.write('%d %d %d\n'%(gain,start_pos_deg,stepsPerRotation))
fd.flush()

dev.set_dio_lo(ledPin) # turn LED off
	
# rotate to start postion
current_pos = dev.get_pos()
dev.soft_ramp_to_pos(current_pos - numpy.mod(current_pos,stepsPerRotation)+start_pos,RAMP_ACCEL)

pstring =  'trial start pos = '+ str(start_pos_deg) + '\n ./' + wavelength + '_' + filterType + '_' + eyeTreatment + '_' + date_time + '.txt \n'
stdscr.addstr(pstring)
stdscr.refresh()
		
time.sleep(PRE_LED_TIME)

#turn LED on
dev.set_dio_hi(ledPin)

time.sleep(PRE_MOTION_TIME)

dev.set_vel_setpt(0)
dev.set_mode('velocity')
dev.start()
start_time = time.time()
t = 0
while t <= start_time + TOTAL_TIME_EXP:
	gin = stdscr.getch()
	if gin ==10:
		stdscr.addstr("experiment ended")
		stdscr.refresh()
		time.sleep(1)
		break
	else:
		rawAnalogIn = comedi.comedi_data_read(daq,subdev,chans[0],0,comedi.AREF_GROUND)
		voltsInLMR = comedi.comedi_to_phys(rawAnalogIn[1], cr[0], maxdata[0])

		rawAnalogIn = comedi.comedi_data_read(daq,subdev,chans[1],0,comedi.AREF_GROUND)
		voltsInFreq = comedi.comedi_to_phys(rawAnalogIn[1], cr[1], maxdata[1])

		vel = int(round(voltsInLMR*gain))

		dev.set_dir_setpt(direction[(numpy.sign(vel)+1)/2], io_update=False)
		dev.set_vel_setpt(abs(vel))

		t = time.time()
		pos = dev.get_pos()

		fd.write('%f %d %d %f %f\n'%(t,vel,pos,voltsInLMR,voltsInFreq))
		fd.flush()

	
dev.stop()
fd.close()

curses.endwin()
	
dev.set_dio_lo(ledPin) # turn LED off
dev.soft_ramp_to_vel(0,'positive',RAMP_ACCEL)
current_pos = dev.get_pos()
dev.soft_ramp_to_pos(current_pos - numpy.mod(current_pos,stepsPerRotation),RAMP_ACCEL)

sys.stdout.write('\a') # make a beep sound to mark end

panels.SetPositions(143, 0)  # load [single] stripe pattern
panels.SetMode(modex, modey)
panels.SetGainOffset(gainx,offsetx,gainy,offsety)
panels.Start()







