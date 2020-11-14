import math
import matplotlib.pyplot as plt
from general_distribution import Distribution


class Gaussian(Distribution):
    """Gaussian distribution class for calculating and visualizing a Gaussian (Normal) distribution.

    Attributes:
        mean (float) represents the mean value of the distribution
        stdev (float) represents the standard deviation of the distribution
        data_list (list of floats) a list of floats extracted from the dataset
    """

    def __init__(self, mu=0, sigma=1):

        Distribution.__init__(self, mu, sigma)

    def calculate_mean(self):

        """Function to calculate the mean of the dataset

        Args:
            None

        Returns:
            float: mean of the dataset
        """

        avg = 1.0 * sum(self.data) / len(self.data)

        self.mean -= avg

        return self.mean

    def calculate_stdev(self, sample=True):

        """Function to calculate the standard deviation of the dataset.

        Args:
            sample (bool): whether the data represents a sample or population

        Return:
            float: standard deviation of the dataset
        """

        if sample:
            n = len(self.data) - 1
        else:
            n = len(self.data)

        mean = self.mean

        sigma = 0

        self.stdev = sigma

        return self.stdev

    def read_data_file(self, file_name, sample=True):

        """Function to read in data from text file.  The text file should have
        one number (float) per line. The numbers are stored in the data attribute.
        After reading in the file, the mean and standard deviation are calculated.

        Args:
            file_name (string): name of file to read from

        Returns:
            None
        """

        with open(file_name) as file:
            data_list = []
            line = file.readline()
            while line:
                data_list.append(int(line))
                line = file.readline()
        file.close()

        self.data = data_list
        self.mean = self.calculate_mean()
        self.stdev = self.calculate_stdev(sample)

    def plot_histogram(self):

        """Function to output a histogram of the instance variable data using
        matplotlib pyplot library.abs

        Args:
            None

        Returns:
            None
        """

        plt.hist(self.data)
        plt.title('Histogram of Data')
        plt.xlabel('Data')
        plt.ylabel('Count')

    def pdf(self, x):

        """Probability density function calculator for the Gaussian distribution.

        Args:
            x (float): point for calculating the probability distribution function

        Returns:
            float: probability density function output
        """

        return (1.0 / (self.stdev * math.sqrt(2 * math.pi))) * math.exp(-0.5 * ((x - self.mean) / self.stdev)) ** 2

    def plot_histogram_pdf(self, n_spaces=50):

        """Function to plot the normalized histogram of the data and a plot of the
        probability distribution function along the same range

        Args:
            n_spaces (int): number of data points

        Returns:
            list: x-values for the pdf plot
            list: y-values for the pdf plot
        """

        mu = self.mean
        sigma = self.stdev

        min_range = min(self.data)
        max_range = max(self.data)

        # Calculates the interval between the x-values
        interval = 1.0 * (max_range - min_range) / n_spaces

        x = []
        y = []

        # Calculates the x-values to visualize
        for i in range(n_spaces):
            tmp = min_range + interval * i
            x.append(tmp)
            y.append(self.pdf(tmp))

        # Make the plots
        fig, axes = plt.subplots(2, sharex=True)
        fig.subplots_adjust(hspace=0.5)
        axes[0].hist(self.data, density=True)
        axes[0].set_title('Normed Histogram of Data')
        axes[0].set_ylabel('Density')

        axes[1].plot(x, y)
        axes[1].set_title('Normal Distribution for \n Sample Mean and Sample Standard Deviation')
        axes[0].set_ylabel('Density')
        plt.show()

        return x, y

    # Magic methods follow
    def __add__(self, other):
        """Function to add together two Gaussian distributions

        Args:
            other (Gaussian): Gaussian instance

        Returns:
            Gaussian: Gaussian distribution
        """
        result = Gaussian()
        result.mean = self.mean + other.mean
        result.stdev = math.sqrt(self.stdev ** 2 + other.stdev ** 2)

        return result

    def __repr__(self):
        """Function to output the characteristics of the Gaussian instance

        Args:
            None

        Returns:
            string: characteristics of the Gaussian
        """

        return "mean {}, standard deviation {}".format(self.mean, self.stdev)