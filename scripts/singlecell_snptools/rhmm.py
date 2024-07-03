import os
import json
from copy import copy
import itertools as it
from collections import Counter

os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np
import pandas as pd
import pomegranate as pm

MAX_ERROR = 0.05


class RigidHMM:

    def __init__(self, rfactor, term_rfactor, trans_prob, error_factor=100):
        self.rfactor = rfactor
        self.term_rfactor = term_rfactor
        self.trans_prob = trans_prob
        self.error_factor = error_factor
        self.error_prob = min(trans_prob * error_factor, MAX_ERROR)
        self.term_prob = 1 / (self.rfactor - self.term_rfactor)
        self._model = None
        self._distributions = None
        self._state_map = None
        self._idx_pos_class = None
        self._initialised = False

    def initialise_model(self, dists):
        self._model = pm.HiddenMarkovModel()
        self._distributions = dists
        start_states = []
        end_states = []
        for s, dist in enumerate(dists):
            prev_state = pm.State(dist, name=f'{s}_0')
            self._model.add_state(prev_state)
            self._model.add_transition(self._model.start, prev_state, self.term_prob)
            self._model.add_transition(prev_state, prev_state, 1 - self.trans_prob)
            start_states.append(prev_state)
            noise_dist = dists[1 - s]
            prev_noise_state = pm.State(noise_dist, name=f'{s}_0_noise')
            self._model.add_transition(prev_state, prev_noise_state, self.error_prob)
            self._model.add_transition(prev_noise_state, prev_state, self.error_prob)
            for i in range(1, self.rfactor):
                state = pm.State(dist, name=f'{s}_{i}')
                self._model.add_state(state)
                self._model.add_transition(prev_state, state, self.trans_prob if i == 1 else 1)
                self._model.add_transition(prev_noise_state, state, 1)
                noise_state = pm.State(noise_dist, name=f'{s}_{i}_noise')
                self._model.add_transition(prev_state, noise_state, self.error_prob)
                # allow shorter states at the termini by adding shortcuts to start/end
                if i <= (self.rfactor - self.term_rfactor):
                    self._model.add_transition(self._model.start, state, self.term_prob)
                if (i + 1) >= self.term_rfactor:
                    self._model.add_transition(state, self._model.end, self.term_prob)
                prev_state = state
                prev_noise_state = noise_state
            else:
                self._model.add_transition(prev_noise_state, prev_state, 1)
                self._model.add_transition(prev_state, prev_noise_state, self.error_prob)
                end_states.append(prev_state)
        for i, j in it.product(range(2), repeat=2):
            if i == j:
                self._model.add_transition(end_states[i], end_states[j], 1 - self.trans_prob)
            else:
                self._model.add_transition(end_states[i], start_states[j], self.trans_prob)
        self._model.bake()
        self._state_map = np.array([
            int(s.name.split('_')[0])
            if not s.name.startswith('None')
            else None
            for s in self._model.states])
        self._idx_pos_class = np.where(self._state_map == 1)[0]
        self._initialised = True

    def fit(self, X, fit_hmm=False, **kwargs):
        init_points = np.concatenate(X)
        assert init_points.shape[1] == 2 and len(init_points.shape) == 2
        fg_poisson = pm.PoissonDistribution.from_samples(init_points.max(1))
        bg_poisson = pm.PoissonDistribution.from_samples(init_points.min(1))
        ref_state = pm.IndependentComponentsDistribution([fg_poisson, bg_poisson])
        alt_state = pm.IndependentComponentsDistribution([bg_poisson, fg_poisson])
        self.initialise_model([ref_state, alt_state])
        if fit_hmm:
            self._model.fit(X, distribution_inertia=0.1, **kwargs)
        else:
            self._model.fit(X, distribution_inertia=0.1, edge_inertia=1, **kwargs)
        return self
    
    def predict(self, markers):
        if not self._initialised:
            raise ValueError('Model has not been initialised')

        # predict in all four orientations and merge
        proba = []
        for hap_flip, seq_flip in it.product(range(2), repeat=2):
            markers_aug = np.flip(markers, axis=1) if hap_flip else markers
            markers_aug = np.flip(markers_aug, axis=0) if seq_flip else markers_aug
            proba_aug = self._model.predict_proba(markers_aug)
            proba_aug = proba_aug[:, self._idx_pos_class].sum(1)
            proba_aug = 1 - proba_aug if hap_flip else proba_aug
            proba_aug = np.flip(proba_aug, axis=0) if seq_flip else proba_aug
            proba.append(proba_aug)
        proba = np.mean(proba, axis=0)
        states = (proba > 0.5).astype(int)
        return states, proba

    def save(self, fn):
        with open(fn, 'w') as o:
            pm_hmm_json = json.loads(self._model.to_json())
            pm_dists_json = [json.loads(d.to_json()) for d in self._distributions]
            full_json = {
                'class': 'RigidHMM',
                'params': {
                    'rfactor': self.rfactor,
                    'term_rfactor': self.term_rfactor,
                    'trans_prob': self.trans_prob,
                    'error_factor': self.error_factor,
                },
                'pm_dists': pm_dists_json,
                'pm_hmm': pm_hmm_json
            }
            json.dump(full_json, o)

    @classmethod
    def load_model(cls, fn):
        with open(fn) as f:
            full_json = json.load(f)
        assert full_json['class'] == 'RigidHMM'
        rhmm = cls(**full_json['params'])
        rhmm._distributions = [pm.NormalDistribution.from_json(json.dumps(d)) for d in full_json['pm_dists']]
        rhmm._model = pm.HiddenMarkovModel.from_json(json.dumps(full_json['pm_hmm']))
        rhmm._state_map = np.array([
            int(s.name.split('_')[0])
            if not s.name.startswith('None')
            else None
            for s in rhmm._model.states])
        rhmm._idx_pos_class = np.where(rhmm._state_map == 1)[0]
        rhmm._initialised = True
        return rhmm


def trim_idx(markers):
    bin_sums = markers.sum(1)
    first = 0
    for i in bin_sums:
        if i != 0.:
            break
        else:
            first = first + 1
    last = len(bin_sums)
    for i in bin_sums[::-1]:
        if i != 0.:
            break
        else:
            last = last - 1
    return first, last


def create_training_data(co_markers, nexamples=100, min_length=50):
    training_data = []
    f_counts = Counter()
    for (f_idx, cb), cb_co_markers in co_markers.items():
        if f_counts[f_idx] >= nexamples:
            continue
        else:
            for x in cb_co_markers.values():
                x = x.astype(np.float32)
                if not x.sum():
                    # empty chrom
                    continue
                s, e = trim_idx(x)
                if (e - s) > min_length:
                    x = x[s:e].copy()
                    if np.random.randint(low=0, high=2):
                        x = np.flip(x, axis=1)
                    if np.random.randint(low=0, high=2):
                        x = np.flip(x, axis=0)
                    training_data.append(x)
                    f_counts[f_idx] += 1
                if f_counts[f_idx] >= nexamples:
                    break
    return training_data