import numpy as np
from scipy.linalg import det
from scipy.linalg import inv


class BhattacharyyaDistance:

    # Distance between multivariate distributions.
    # group_one, group_two are lists that contain
    # composite elements of the same dimension
    # e.g. for 2 dimensions: [ [1,2], [2,3], [6.1] ]
    @staticmethod
    def multivariate_compare(group_one, group_two):
        group_one_mean = np.mean(group_one, axis=0)
        group_two_mean = np.mean(group_two, axis=0)
        group_one_cov = np.cov(group_one.T)
        group_two_cov = np.cov(group_two.T)
        aggregate_cov = (group_one_cov + group_two_cov) / 2
        mean_difference = group_one_mean - group_two_mean
        mahalanobis_term = np.dot(np.dot(mean_difference.T, inv(aggregate_cov)), mean_difference)
        det_term = det(aggregate_cov) / np.sqrt(det(group_one_cov) * det(group_two_cov))

        return (0.125 * mahalanobis_term) + (0.5 * np.log(det_term))

    @staticmethod
    def multivariate_compare_group(group_one, group_two):
        if(group_one['det_cov'] <= 0 or group_two['det_cov'] <= 0):
            return False
        aggregate_cov = (group_one['cov'] + group_two['cov']) / 2
        mean_difference = group_one['mean'] - group_two['mean']
        mahalanobis_term = np.dot(np.dot(mean_difference.T, inv(aggregate_cov)), mean_difference)
        det_term = det(aggregate_cov) / np.sqrt(group_one['det_cov'] * group_two['det_cov'])

        return (0.125 * mahalanobis_term) + (0.5 * np.log(det_term))
