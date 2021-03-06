**Enabling the watchdog service**

Sometimes the RPi whill hang for unknown reasons, so it is a good idea to restart it then. This can be done automatically using the onboard watchdog timer.

1. Install the watchdog package:

```
sudo apt-get install watchdog
```

2. Load the module manually:

```
sudo modprobe bcm2835_wdt
```

3. Then, add a config to automatically load the module:
Open config file:
```
sudo nano /etc/modules-load.d/bcm2835_wdt.conf
```

Add this line in the file and save it:

```
bcm2835_wdt
```

4. We’ll also need to manually edit the systemd unit at /lib/systemd/system/watchdog.service and add a line to the [Install] section:

```
[Install]
WantedBy=multi-user.target
```

5. We need to configure the watchdog.

Open /etc/watchdog.conf with your favorite editor.

```
sudo nano /etc/watchdog.conf
```

Uncomment the line that starts with #watchdog-device by removing the hash (#) to enable the watchdog daemon to use the watchdog device.

Uncomment the line that says #max-load-1 = 24 by removing the hash symbol to reboot the device if the load goes over 24 over 1 minute. A load of 25 of one minute means that you would have needed 25 Raspberry Pis to complete that task in 1 minute. You may tweak this value to your liking.

6. Then enable the service:

```
sudo systemctl enable watchdog.service
```

7. Finally, start the service:

```
sudo systemctl start watchdog.service
```

8. You can set various options for the watchdog in /etc/watchdog.conf – see the man page for that file.
