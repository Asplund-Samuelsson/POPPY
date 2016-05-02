# Progress indicator for Python 3.5.1
# Version 0.1.1
# Johannes Asplund-Samuelsson

class Progress():
    """Progress indication thingamabob"""

    # Design variables
    __allowed_designs = set(['p','s','b','t','c'])

    # Spinner variables
    __spin_stages = ['/','-','\\','|']
    __spin_stage = 0
    __last_spin_time = 0

    # Timer variables
    __previous_progress = 0
    __last_timer_time = 0

    def __init__(self, max_val = 100, design = 'p', val = 0):
        """Initialize the Progress indicator

        ARGUMENTS

        max_val : int, float
            The expected maximum or target value.
        design : string
            The type and order of progress indicators.
            'p' : percent
            's' : spinner
            'b' : bar
            't' : timer
            'c' : counter
        val : int, float
            The current value.
        """

        # Import functions
        from os import popen
        from sys import stdout
        from time import time
        from time import sleep
        from collections import deque
        from datetime import timedelta

        # Bind imported functions to self
        self.__popen = popen
        self.__stdout = stdout
        self.__time = time
        self.__sleep = sleep
        self.__deque = deque
        self.__timedelta = timedelta

        # Initialize variables
        self.update(val, max_val, design)
        self.__speed_samples = self.__deque(maxlen=10)

    def __call__(self, val = None):
        """Calling returns a string of the current progress"""
        if val:
            self.val = val
        return self.to_string()

    def __format(self):
        """Format the progress indicator output string"""
        output = []
        for variant in self.design:
            if variant == 'p':
                output.append(self.percent())
            if variant == 's':
                output.append(self.spinner())
            if variant == 'b':
                output.append(self.bar())
            if variant == 't':
                output.append(self.timer())
            if variant == 'c':
                output.append(self.counter())
        return ' '.join(output)

    def __avg_speed(self):
        """Calculate the average speed"""
        if self.__speed_samples:
            avg_speed = sum([s[0] / s[1] for s in self.__speed_samples]) \
            / len(self.__speed_samples)
        else:
            avg_speed = 0
        return avg_speed

    def percent(self):
        """Percent progress towards the maximum"""
        # Calculate percent as float and format as string
        if self.max_val:
            percent = float(self.val / self.max_val * 100)
            return "{0:>6}".format("%0.1f%%" % percent)
        else:
            percent = "INF%"
            return "{0:>6}".format(percent)

    def spinner(self):
        """A spinner that (might) indicate that something is happening"""
        # Calculate time since last spin
        time_since_spin = self.__time() - self.__last_spin_time
        # Determine whether to spin or not
        # Don't spin...
        if time_since_spin < 0.1:
            return self.__spin_stages[self.__spin_stage]
        # Spin...
        # Calculate spin stage
        if self.__spin_stage < len(self.__spin_stages) - 1:
            self.__spin_stage += 1
        else:
            self.__spin_stage = 0
        # Set time of last spin
        self.__last_spin_time = self.__time()
        # Return spinner character of current spin stage
        return self.__spin_stages[self.__spin_stage]

    def bar(self):
        """Progress bar with the wget design"""
        # Get the current terminal width
        try:
            rows, columns = self.__popen('stty size', 'r').read().split()
            columns = int(columns)
        except ValueError:
            columns = 80

        # Count the width of other indicators
        width = 0
        for variant in self.design:
            if variant == 'p':
                width += len(self.percent())
                width += 1
            if variant == 's':
                width += 2
            if variant == 't':
                width += len(self.timer())
                width += 1
            if variant == 'c':
                width += len(self.counter())

        # Calculate the allowed bar width
        allowed_width = min(columns - width, 55)

        # Calculate the filled width at current progress state
        full_width = allowed_width - 2 # 2 for the caps

        # If the maximum value is zero, the bar cannot be calculated
        if not self.max_val:
            return ''.join(['['] + ['='] * full_width + [']'])

        # Construct the progress bar - Example: [=====>       ]
        n_filled = round((self.val / self.max_val) * full_width)
        n_empty = full_width - n_filled
        bar = ['='] * n_filled
        if bar and n_empty > 0:
            bar[-1] = '>'
        bar = ['['] + bar + [' '] * n_empty + [']']

        return ''.join(bar)

    def timer(self):
        """Countdown timer based on average progress speed"""
        time_since_timer = self.__time() - self.__last_timer_time
        # Add a speed sample if 1 or more seconds have passed
        if time_since_timer >= 1:
            self.__speed_samples.append(
                (self.val - self.__previous_progress, time_since_timer)
            )
            self.__previous_progress = self.val
            self.__last_timer_time = self.__time()
        # If an average speed can be calculated, return a remaining time string
        if self.__avg_speed():
            remaining_time = (self.max_val - self.val) / self.__avg_speed()
            if remaining_time <= 24*3600-1:
                clock = str(self.__timedelta(seconds=int(remaining_time)))
                return "{0:>8}".format(clock)
            elif remaining_time <= 7*24*3600-1:
                days = str(self.__timedelta(seconds=int(remaining_time)).days)
                return ' >' + days + ' days'
            else:
                return ' >1 week'
        # Return a blank remaining time string if the average speed is zero
        else:
            return '--:--:--'

    def counter(self):
        """The actual value as a fraction of the maximum"""
        form = "{0:>%s}" % str(2*len(str(self.max_val)) + 2)
        return form.format(str(self.val) + '/' + str(self.max_val))

    def update(self, val = None, max_val = None, design = None):
        """Update parameters of the progress indicator"""
        if val is not None:
            self.val = val
        if max_val is not None:
            self.max_val = max_val
        if design is not None:
            if set(design) - self.__allowed_designs:
                raise Exception('Invalid progress indicator type.')
            else:
                self.design = list(design)

    def to_string(self, val = None):
        """Produce a string of the current progress"""
        if val:
            self.val = val
        return self.__format()

    def write(self, val = None):
        """Write progress to standard output with carriage return and flush."""
        if val:
            self.val = val
        self.__stdout.write("\r" + self.__format())
        self.__stdout.flush()

    def test(self):
        """Test the output of the active parameters"""
        val = self.val
        for i in range(1000):
            x = (i + 1) / 10
            self.write(x)
            self.__sleep(0.01)
        self.val = val
        self.__stdout.write("\n")
