# RPi Meteor Station
# Copyright (C) 2015  Dario Zubovic
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from multiprocessing import Process, Event
import cv2
import time
import logging

# Get the logger from the main module
log = logging.getLogger("logger")


class BufferedCapture(Process):
    """ Capture from device to buffer in memory.
    """
    
    running = False
    
    def __init__(self, array1, timearray1, startTime1, array2, timearray2, startTime2, config, video_file=None):
        """ Populate arrays with (startTime, frames) after startCapture is called.
        
        Arguments:
            array1: numpy array in shared memory that is going to be filled with frames
            timearray1: numpy array in shared memory containing time stamps
            startTime1: float in shared memory that holds time of first frame in array1
            array2: second numpy array in shared memory
            timearray2: second numpy array in shared memory containing time stamps
            startTime2: float in shared memory that holds time of first frame in array2

        Keyword arguments:
            video_file: [str] Path to the video file, if it was given as the video source. None by default.

        """
        
        super(BufferedCapture, self).__init__()
        self.array1 = array1
        self.timearray1 = timearray1
        self.startTime1 = startTime1
        self.array2 = array2
        self.timearray2 = timearray2
        self.startTime2 = startTime2
        
        self.startTime1.value = 0
        self.startTime2.value = 0
        
        self.config = config

        self.video_file = video_file

        # A fraem will be considered dropped if it was late more then half a frame
        self.time_for_drop = 1.5*(1.0/config.fps)

        self.dropped_frames = 0
    


    def startCapture(self, cameraID=0):
        """ Start capture using specified camera.
        
        Arguments:
            cameraID: ID of video capturing device (ie. ID for /dev/video3 is 3). Default is 0.
            
        """
        
        self.cameraID = cameraID
        self.exit = Event()
        self.start()
    


    def stopCapture(self):
        """ Stop capture.
        """
        
        self.exit.set()
        self.join()



    def run(self):
        """ Capture frames.
        """
        
        # use a file as the video source
        if self.video_file is not None:
            device = cv2.VideoCapture(self.video_file)

        # Use a device as the video source
        else:

            # Init the video device
            device = cv2.VideoCapture(self.config.deviceID)

            # Try setting the resultion if using a video device, not gstreamer
            try:

                # This will fail if the video device is a gstreamer pipe
                int(self.config.deviceID)
                
                # Set the resolution (crashes if using an IP camera and gstreamer!)
                device.set(3, self.config.width)
                device.set(4, self.config.height)

            except:
                pass


        # Wait until the device is opened
        device_opened = False
        for i in range(20):
            time.sleep(0.1)
            if device.isOpened():
                device_opened = True
                break


        # If the device could not be opened, stop capturing
        if not device_opened:
            log.info('The video source could not be opened!')

            return False

        else:
            log.info('Video device opened!')



        # Throw away first 10 frame
        for i in range(10):
            device.read()


        first = True
        
        # Run until stopped from the outside
        while not self.exit.is_set():
            lastTime = 0
            
            if first:
                self.startTime1.value = 0
            else:
                self.startTime2.value = 0
            
            for i in range(256):

                # Read the frame
                ret, frame = device.read()

                # If the end of the file was reached, stop the capture
                if (self.video_file is not None) and (frame is None):

                    log.info('End of video file! Press Ctrl+C to finish.')

                    self.exit.set()
                    
                    time.sleep(0.1)

                    break

                
                t = time.time()

                if i == 0: 
                    startTime = t

                # check if frame is dropped
                elif lastTime - t >= self.time_for_drop: 
                    
                    # Calculate the number of dropped frames
                    n_dropped = int((lastTime - t)*self.config.fps)
                    
                    log.info(str(n_dropped) + " frames dropped!")

                    self.dropped_frames += n_dropped

                    

                lastTime = t
                
                # Convert the frame to grayscale
                gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

                if first:
                    self.array1[i, :gray.shape[0], :gray.shape[1]] = gray
                    self.timearray1[i] = t
                else:
                    self.array2[i, :gray.shape[0], :gray.shape[1]] = gray
                    self.timearray2[i] = t


                # If video is loaded from a file, simulate real FPS
                if self.video_file is not None:

                    time.sleep(1.0/self.config.fps)

                    # If the video finished, stop the capture
                    if not device.isOpened():
                        self.exit.set()



            if self.exit.is_set():
                break


            if first:
                self.startTime1.value = startTime

            else:
                self.startTime2.value = startTime

            first = not first
        

        device.release()
    
