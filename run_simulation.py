import numpy as np
import matplotlib.pyplot as plt

from pylocus.basics import vector_from_matrix
from cvxpy import *

#FOLDER = 'gaussian_results_corrected'
FOLDER = 'reproduce_test'


def get_noisy_inner(Om_original, dm_original, rho, sigmad):
    from tikhonov import tikhonov_scipy
    tik_sc = tikhonov_scipy(rho)
    errors = np.triu(tik_sc.rvs(size=Om_original.shape))
    Om = np.cos(np.arccos(Om_original) + errors + errors.T)
    np.fill_diagonal(Om, 1.0)
    dm = (np.sqrt(dm_original) + np.random.normal(0,
                                                  sigmad, Om_original.shape[0]))**2
    outer = np.outer(dm, dm)
    KE = np.multiply(outer, Om)
    return Om, dm, KE


def get_Om_from_abs_angles(abs_angles, m):
    angles_vector = vector_from_matrix(abs_angles)
    Om = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            Om[i, j] = np.cos(angles_vector[i] - angles_vector[j])
    return Om


def get_noisy_absolute(abs_angles_original, dm_original, rho, sigmad, gaussian):
    m = len(dm_original)
    if gaussian:
        errors = np.triu(np.random.normal(0, rho, abs_angles_original.shape))
    else:
        from tikhonov import tikhonov_scipy
        tik_sc = tikhonov_scipy(rho)
        errors = np.triu(tik_sc.rvs(size=abs_angles_original.shape))
    abs_angles_noisy = abs_angles_original + errors + errors.T
    np.fill_diagonal(abs_angles_noisy, 0.0)
    abs_angles_vector = vector_from_matrix(abs_angles_noisy)

    if gaussian:
        # create Om with same variance as absolute angles.
        errors = np.triu(np.random.normal(
            0, rho/2.0, abs_angles_original.shape))
        abs_angles_noisy = abs_angles_original + errors + errors.T
    Om = get_Om_from_abs_angles(abs_angles_noisy, m)
    dm = (np.sqrt(dm_original) + np.random.normal(0, sigmad, m))**2
    return Om, dm, abs_angles_vector


def get_noisy_points(points_original, noise):
    from pylocus.plots_cti import plot_matrix
    from pylocus.basics_angles import from_0_to_2pi
    points_noisy = points_original.copy()
    points_noisy.add_noise(noise)
    abs_angles_vector = vector_from_matrix(points_noisy.abs_angles)
    #diff_angles = points_noisy.abs_angles - points_original.abs_angles
    #diff_dm = (points_noisy.dm - points_original.dm).reshape((1,-1))
    #diff = from_0_to_2pi(diff)
    #plot_matrix(diff_angles, 'diff abs angles')
    #print('diff dm:',diff_dm)
    return points_noisy.Om, points_noisy.dm, abs_angles_vector


def parse_options(options):
    from pylocus.point_set import HeterogenousSet
    sigmas = np.linspace(options['min_sigma'],
                         options['max_sigma'], options['n_sigma'])
    rhos = np.exp(np.linspace(
        options['min_rho'], options['max_rho'], options['n_rhos']))
    points = HeterogenousSet(options['N'], options['d'])
    return sigmas, rhos, points


def parse_options_gaussian(options):
    from pylocus.point_set import HeterogenousSet
    sigmas = np.linspace(options['min_sigma'],
                         options['max_sigma'], options['n_sigma'])
    rhos = np.linspace(options['min_rho'],
                       options['max_rho'], options['n_rhos'])
    points = HeterogenousSet(options['N'], options['d'])
    return sigmas, rhos, points


def clean_angles(Om, print_out):
    E = len(Om)
    X = Semidef(E)
    Noise = Variable(E, E)

    constraints = [X + Noise == Om]  # include weighting matrix?
    [constraints.append(X[i, i] == 1.0) for i in range(E)]

    obj = Minimize(trace(X) + norm(Noise))
    prob = Problem(obj, constraints)

    #total_cost = prob.solve(solver='SCS',verbose=True, eps=1e-10)
    #total_cost = prob.solve(solver='CVXOPT',verbose=True,
    #                        abstol=1e-10, reltol=1e-8, feastol=1e-10,
    #                        kktsolver="robust")
    total_cost = prob.solve(solver='CVXOPT', verbose=print_out,
                            kktsolver="robust")
    if X.value is not None:
        return X.value, Noise.value


def run_simulation(methods, options, save_idx=None):
    from pylocus.algorithms import reconstruct_emds, reconstruct_mds, reconstruct_smds
    from pylocus.basics import rmse
    from pylocus.point_set import edm_from_dm
    if options['gaussian']:
        sigmas, rhos, points = parse_options_gaussian(options)
    else:
        sigmas, rhos, points = parse_options(options)
    print_out = options['print_out']

    points.set_points('normal')

    C, b = points.get_KE_constraints()

    dict_methods = {m: {'rmses': np.zeros(
        (len(sigmas), len(rhos))), 'estimate': ''} for m in methods}

    for n in range(options['n_it']):
        points.set_points('normal')
        print('n', n)
        for j, rho in enumerate(rhos):
            print('    angle noise j', j)
            for i, sigmad in enumerate(sigmas):
                print('      distance noise i', i)
                Om, dm, absolute_angles = get_noisy_absolute(
                    points.abs_angles, points.dm, rho, sigmad, options['gaussian'])
                edm = edm_from_dm(dm, points.N)
                # TODO only treat chosen methods.
                for key in dict_methods.keys():
                    if key == 'E-MDS':
                        dict_methods[key]['estimate'] = reconstruct_emds(
                            edm, real_points=points.points, Om=Om)
                    elif key == 'MDS':
                        dict_methods[key]['estimate'] = reconstruct_mds(
                            edm, real_points=points.points)
                    elif key == 'constrained E-MDS':
                        dict_methods[key]['estimate'] = reconstruct_emds(
                            edm, real_points=points.points, Om=Om, method='iterative', C=C, b=b)
                    elif key == 'relaxed E-MDS':
                        dict_methods[key]['estimate'] = reconstruct_emds(
                            edm, real_points=points.points, Om=Om, method='relaxed', C=C, b=b)
                    elif key == 'enhanced E-MDS':
                        clean_Om, Noise = clean_angles(Om, print_out)
                        if print_out:
                            print('difference clean_Om, Om',
                                  np.linalg.norm(clean_Om - Om))
                        dict_methods[key]['estimate'] = reconstruct_emds(
                            edm, real_points=points.points, Om=clean_Om)
                    elif key == 'CDM':
                        dict_methods[key]['estimate'] = reconstruct_smds(
                            dm, absolute_angles, real_points=points.points)
                    else:
                        raise NameError('Unknown method', key)
                    dict_methods[key]['rmses'][i, j] += rmse(
                        dict_methods[key]['estimate'], points.points)
                    if print_out:
                        print('rmse {}: {}'.format(
                            key, dict_methods[key]['rmses'][i, j]))

    if (save_idx):
        import json
        with open('{}/options_{}.json'.format(FOLDER, save_idx), 'w') as outfile:
            json.dump(options, outfile)
        print('saved {}/options_{}'.format(FOLDER, save_idx))
        for m in dict_methods.keys():
            dict_methods[m]['rmses'] /= options['n_it']
            np.save('{}/rmses_{}_{}'.format(FOLDER, m, save_idx),
                    dict_methods[m]['rmses'])
            print('saved {}/rmses_{}_{}'.format(FOLDER, m, save_idx))


if __name__ == '__main__':
    import time

    # Choose simulation paramters.
    options = {
        'N': 5,  # number of points to localize
        'd': 2,  # dimension of points
        'n_sigma': 3,  # was 10, number of distance noise levels to run for
        'n_it': 2, # number of point sets to average over. (100 for smooth curves)
        'n_rhos': 4,  # number of angle noises to test for.
        'min_sigma': 0.01,  # minimum distance noise
        'max_sigma': 0.5,  # maximum distance noise
        'min_rho': 0.01,  # minimum angle noise
        'max_rho': 0.5,  # maximum angle noise
        'print_out': False,  # print debugging information
        'gaussian': True  # if Gaussian noise or Tikhonov noise should be used for angles
    }

    # Choose methods to run.
    # can be any subset of:
    #    - CDM (coordinate difference matrices, described in paper)
    #    - E-MDS (edge-kernel method)
    #    - constrainted E-MDS (constrained edge-kernel method)
    #    - MDS (multidimensional scaling, uses distances only)
    methods = ['CDM', 'E-MDS', 'constrained E-MDS', 'MDS']  # all
    #methods = ['relaxed E-MDS']

    # Run simulations.
    save_idx = int(time.time())  # appendix to add to name that is being saved
    #save_idx = None
    run_simulation(methods, options, save_idx)
