import json
from copy import copy
import itertools as it
import numpy as np
import pandas as pd
import pomegranate as pm


class RigidHMM:

    def __init__(self, rfactor, term_rfactor, nstates=2, noise=1e-2, filter_tiny=1_000_000):
        self.rfactor = rfactor
        self.term_rfactor = term_rfactor
        self.nstates = nstates
        self.noise = noise
        self.filter_tiny = filter_tiny
        self._model = None
        self._distributions = None
        self._state_map = None
        self._initialised = False

    def initialise_model(self, dists, freeze_dists=True):
        self._model = pm.HiddenMarkovModel()
        self._distributions = dists
        start_states = []
        end_states = []
        for s, dist in enumerate(dists):
            if freeze_dists:
                dist.frozen = True
            prev_state = pm.State(dist, name=f'{s}_0')
            self._model.add_state(prev_state)
            self._model.add_transition(self._model.start, prev_state, 1 / self.nstates)
            start_states.append(prev_state)
            for i in range(1, self.rfactor):
                state = pm.State(dist, name=f'{s}_{i}')
                self._model.add_state(state)
                self._model.add_transition(prev_state, state, 1)
                # allow shorter states at the termini by adding shortcuts to start/end
                if i <= (self.rfactor - self.term_rfactor):
                    self._model.add_transition(self._model.start, state, 1 / self.nstates)
                if (i + 1) >= self.term_rfactor:
                    self._model.add_transition(state, self._model.end, 1 / self.nstates)
                prev_state = state
            end_states.append(prev_state)
        for i, j in it.product(range(self.nstates), repeat=2):
            if i == j:
                tprob = 0.9
                self._model.add_transition(end_states[i], end_states[j], tprob)
            else:
                tprob = 0.1 / (self.nstates - 1)
                self._model.add_transition(end_states[i], start_states[j], tprob)
        self._model.bake()
        self._state_map = np.array([
            int(s.name.split('_')[0])
            if not s.name.startswith('None')
            else None
            for s in self._model.states])
        self._initialised = True

    def fit(self, X, **kwargs):
        gmm = pm.GeneralMixtureModel.from_samples(pm.NormalDistribution, self.nstates, np.concatenate(X))
        dists = sorted(gmm.distributions, key=lambda dist: dist.parameters[0])
        self.initialise_model(dists)
        self._model.fit(X, distribution_inertia=1.0, **kwargs)
        return self

    def _ml_adjust_boundaries(self, sequence, states, positions):
        chgp = np.where(np.diff(states))[0] + 1
        chgp = np.insert(chgp, len(chgp), len(states))
        chgp_pos = []
        chgp_left_marker = []
        chgp_right_marker = []
        left_state = states[0]
        chgp_states = [states[0]]
        left_i = 0
        for i in range(len(chgp) - 1):
            right_i = chgp[i + 1]
            right_state = 1 - left_state
            best_score = -np.inf
            best_j = None
            for j in range(left_i + 1, right_i):
                score = self._distributions[left_state].probability(sequence[left_i: j]).sum() + \
                        self._distributions[right_state].probability(sequence[j: right_i]).sum()
                if score > best_score:
                    best_j = j
                    best_score = score
            left_i = best_j
            chgp_states.append(right_state)
            left_state = right_state
            chgp_pos.append((positions[best_j - 1] + positions[best_j]) // 2)
            chgp_left_marker.append(positions[best_j - 1])
            chgp_right_marker.append(positions[best_j])
        return chgp_pos, chgp_left_marker, chgp_right_marker, chgp_states


    def _eval_sequence_prediction(self, pos, markers, seg_starts, seg_ends, states):
        segs = zip(seg_starts, seg_ends, states)
        start, end, state = next(segs)
        logprobs = []
        nmarkers = []
        state_logprob = 0
        state_nmarkers = 0
        for p, m in zip(pos, markers):
            while True:
                if start < p <= end:
                    state_logprob += self._distributions[state].log_probability(m)
                    state_nmarkers += 1
                    break
                elif p > end:
                    logprobs.append(state_logprob)
                    nmarkers.append(state_nmarkers)
                    state_nmarkers = 0
                    state_logprob = 0
                    start, end, state = next(segs)
                elif p < start:
                    raise ValueError('pos are not sorted')
        else:
            logprobs.append(state_logprob)
            nmarkers.append(state_nmarkers)
        return logprobs, nmarkers


    def _filter_tiny_segments(self, seg_pos, seg_left_marker, seg_right_marker, states, filter_tiny):
        # filter tiny COs
        keep_nontiny = np.diff(seg_pos) > filter_tiny
        # allow chrom ends to be tiny
        keep_nontiny[0] = True
        keep_nontiny[-1] = True
        states = states[keep_nontiny]
        keep_nontiny = np.insert(keep_nontiny, len(keep_nontiny), True)
        seg_pos = seg_pos[keep_nontiny]
        seg_left_marker = seg_left_marker[keep_nontiny]
        seg_right_marker = seg_right_marker[keep_nontiny]
        # now merge adjacent segments with the same state
        keep_diffstate = np.diff(states) != 0
        keep_diffstate = np.insert(keep_diffstate, 0, True)
        states = states[keep_diffstate]
        keep_diffstate = np.insert(keep_diffstate, len(keep_diffstate), True)
        seg_pos = seg_pos[keep_diffstate]
        seg_left_marker = seg_left_marker[keep_diffstate]
        seg_right_marker = seg_right_marker[keep_diffstate]
        return seg_pos, seg_left_marker, seg_right_marker, states
    
    
    def predict(self, pos, markers, chrom_size):
        if not self._initialised:
            raise ValueError('Model has not been initialised')
        states = self._state_map[self._model.predict(markers)]
        seg_pos, seg_left_marker, seg_right_marker, states = self._ml_adjust_boundaries(markers, states, pos)
        seg_pos = np.insert(seg_pos, [0, len(seg_pos)], [0, chrom_size]).astype(np.int32)
        seg_left_marker = np.insert(seg_left_marker, [0, len(seg_left_marker)], [0, chrom_size])
        seg_right_marker = np.insert(seg_right_marker, [0, len(seg_right_marker)], [0, chrom_size])
        states = np.array(states, dtype=np.int32)
        seg_pos, seg_left_marker, seg_right_marker, states = self._filter_tiny_segments(
            seg_pos, seg_left_marker, seg_right_marker, states, self.filter_tiny
        )
        seg_starts = seg_pos[:-1]
        seg_start_min = seg_left_marker[:-1]
        seg_ends = seg_pos[1:]
        seg_end_max = seg_right_marker[1:]
        logprobs, nmarkers = self._eval_sequence_prediction(pos, markers, seg_starts, seg_ends, states)
        return seg_starts, seg_ends, seg_start_min, seg_end_max, states, logprobs, nmarkers

    def save(self, fn):
        with open(fn, 'w') as o:
            pm_hmm_json = json.loads(self._model.to_json())
            pm_dists_json = [json.loads(d.to_json()) for d in self._distributions]
            full_json = {
                'class': 'RigidHMM',
                'params': {
                    'rfactor': self.rfactor,
                    'term_rfactor': self.term_rfactor,
                    'nstates': self.nstates,
                    'noise': self.noise,
                    'filter_tiny': self.filter_tiny,
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
        rhmm._initialised = True
        return rhmm


def get_haplotype_blocks(co_invs):
    haplotype_blocks = {}
    chrom_sizes = co_invs.groupby('chrom').end.max().to_dict()

    for chrom, chrom_co_invs in co_invs.groupby('chrom'):
        haplo_state = (chrom_co_invs.query('start == 0')
                                    .set_index('cb')                      
                                    .haplo.to_dict())
        prev_pos = 0
        chrom_co_invs = chrom_co_invs.query('start > 0')
        for pos, invs in chrom_co_invs.groupby('start', sort=True):
            haplotype_blocks[(chrom, prev_pos, pos)] = copy(haplo_state)
            prev_pos = pos
            # mutate haplo state
            for _, i in invs.iterrows():
                haplo_state[i.cb] = i.haplo
        else:
            haplotype_blocks[(chrom, prev_pos, chrom_sizes[chrom])] = copy(haplo_state)
    return pd.DataFrame.from_dict(haplotype_blocks, orient='index')