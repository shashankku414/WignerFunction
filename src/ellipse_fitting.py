import numpy as np
from scipy.optimize import minimize


class ellipse_fit:
    def __init__(self, Q_values, P_values, W_values):
        self.Q_values = Q_values
        self.P_values = P_values
        self.W_values = W_values
        # self.timedelay = timedelay
    
    def show(self):
        print(self.W_values)

    # this function find the maximum value and its coordinates
    def find_center(self):
        x, y = np.meshgrid(self.Q_values, self.P_values)
        coodinates = list(zip(x.ravel(), y.ravel(), self.W_values.ravel()))

        return max(coodinates, key=lambda t: t[2])

    # this function finds the intersection points between a plane (at value of 1/e * maximum)
    # and the wigner plot
    def find_intersection(self):
        x, y = np.meshgrid(self.Q_values, self.P_values)
        z = self.W_values

        _, _, w_max = self.find_center()
        a, b, c = 0.0, 0.0, w_max/np.exp(1)
        z_plane = a*x + b*y + c

        epsilon = 0.05
        intersection_mask = np.abs(z - z_plane)/z_plane < epsilon

        return x[intersection_mask], y[intersection_mask], z[intersection_mask], z_plane
    
    # def cost_function(self, params, x_data, y_data):
    #     A,B,C,D,E,F = params
    #     err = A*x_data**2 + B*x_data*y_data + C*y_data**2 + D*x_data + E*y_data + F
    #     return np.sum(err**2)
    def cost_function(self, params, x_data, y_data):
        a,b,x0,y0,tilt = params
        err = ((x_data-x0)*np.cos(tilt) + (y_data-y0)*np.sin(tilt))**2/a**2 + ((x_data-x0)*np.sin(tilt) - (y_data-y0)*np.cos(tilt))**2/b**2 - 1
        return np.sum(err**2)
    
    def find_opt_params(self, initial_guess):

        x_intersect, y_intersect, _, _ = self.find_intersection()
        result = minimize(self.cost_function, initial_guess, args=(x_intersect, y_intersect), method='TNC')

        return result.x






# if __name__ == "__main__":
#     obj = ellipse_fit([1,2,3], [2,3,4], [1,1,1])
#     obj.show()