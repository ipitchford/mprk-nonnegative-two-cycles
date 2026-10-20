#!/usr/bin/env python3
"""Exact certificate of an attracting MPRK22(4/7) two-cycle.

Python standard library only. All decisions use rational or integer arithmetic.
The autonomous ODE has A=[[-9990,10],[9990,-10]], h=1 and unit mass.
With t=y_1, u=(7*t+40)/40007, r=t/u, Q(r)=(r+7)*r**(3/4)/8,
the scalar time step is F(t)=(t+10*Q((1-t)/(1-u)))/
(1+9990*Q(t/u)+10*Q((1-t)/(1-u))).

The certificate proves that two disjoint rational intervals map strictly
into one another and that F^2 contracts by a factor below 24/25. Thus the
rational initial first component 1/2 converges to an attracting prime-period
two orbit. A secondary intermediate-value check proves existence independently
of the derivative bounds. No orbit simulation or floating point comparison
enters either proof.
"""
from fractions import Fraction as Q
from math import isqrt
from dataclasses import dataclass
import json
from pathlib import Path

SCALE = 10**30
ROOT_CHECKS = 0

def require(condition):
    if not condition:
        raise ArithmeticError("Exact certificate condition failed")

@dataclass(frozen=True)
class I:
    lo: Q
    hi: Q
    def __init__(self, lo, hi=None):
        object.__setattr__(self, 'lo', Q(lo))
        object.__setattr__(self, 'hi', Q(lo if hi is None else hi))
        require(self.lo <= self.hi)
    def __add__(self, other):
        other = as_i(other)
        return I(self.lo+other.lo, self.hi+other.hi)
    __radd__ = __add__
    def __neg__(self): return I(-self.hi, -self.lo)
    def __sub__(self, other): return self + (-as_i(other))
    def __rsub__(self, other): return as_i(other) + (-self)
    def __mul__(self, other):
        other = as_i(other)
        z=[self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi]
        return I(min(z),max(z))
    __rmul__ = __mul__
    def __truediv__(self, other):
        other = as_i(other)
        require(not (other.lo <= 0 <= other.hi))
        return self*I(1/other.hi,1/other.lo)
    def __rtruediv__(self, other): return as_i(other)/self
    def fourth_root(self):
        require(self.lo >= 0)
        def lower(x):
            global ROOT_CHECKS
            k=isqrt(isqrt((x.numerator*SCALE**4)//x.denominator))
            lo=Q(k,SCALE);hi=Q(k+1,SCALE)
            require(lo**4 <= x < hi**4)
            ROOT_CHECKS += 1
            return lo,hi
        return I(lower(self.lo)[0],lower(self.hi)[1])
    def record(self):
        return {'lower':str(self.lo),'upper':str(self.hi),
                'decimal_lower_display_only':float(self.lo),
                'decimal_upper_display_only':float(self.hi)}

def as_i(x): return x if isinstance(x,I) else I(x)
def weight(r): return (r+7)*(r*r*r).fourth_root()/8

def F(t):
    t=as_i(t)
    require(0 < t.lo <= t.hi < 1)
    u=(7*t+40)/40007
    q1=weight(t/u)
    q2=weight((1-t)/(1-u))
    return (t+10*q2)/(1+9990*q1+10*q2)

def DF(t):
    """Enclose F'(t); differentiate the exact scalar map analytically."""
    t=as_i(t)
    require(0 < t.lo <= t.hi < 1)
    u=(7*t+40)/40007
    r1=t/u;r2=(1-t)/(1-u)
    dw=lambda r: 7*(r+3)/(32*r.fourth_root())
    q1=weight(r1);q2=weight(r2)
    dq1=dw(r1)*Q(40,40007)/(u*u)
    dq2=-dw(r2)*Q(39960,40007)/((1-u)*(1-u))
    num=t+10*q2;den=1+9990*q1+10*q2
    return (1+10*dq2)/den-num*(9990*dq1+10*dq2)/(den*den)

def main():
    global ROOT_CHECKS
    ROOT_CHECKS=0
    left=I(Q(1,4));right=I(Q(3,4));bracket=I(Q(1,4),Q(3,4))
    left_step=F(left);right_step=F(right)
    left_second=F(left_step);right_second=F(right_step)
    whole_step=F(bracket)
    require(left_second.lo > Q(1,4))
    require(right_second.hi < Q(3,4))
    require(whole_step.hi < Q(1,1000))
    require(F(I(Q(1,1000))).lo <= Q(1,1000) <= F(I(Q(1,1000))).hi)
    # Exact local scalar derivative from the published stability function.
    alpha=Q(4,7);z=Q(-10000)
    local_R=(-z*z-2*alpha*z+2)/(2*(1-alpha*z)*(1-z))
    require(abs(local_R) < 1)
    # A concrete rational starting point and an open basin: the two intervals
    # map strictly into one another, and F^2 is a contraction on either one.
    low=I(Q(1,10**7),Q(117,10**9))
    high=I(Q(12,25),Q(507,1000))
    image_low=F(low);image_high=F(high)
    require(high.lo < image_low.lo <= image_low.hi < high.hi)
    require(low.lo < image_high.lo <= image_high.hi < low.hi)
    require(high.lo < Q(1,2) < high.hi)
    derivative_low=DF(low);derivative_high=DF(high)
    require(derivative_low.hi < 0 and derivative_high.hi < 0)
    contraction_bound=derivative_low.lo*derivative_high.lo
    require(0 < contraction_bound < Q(24,25) < 1)
    # These outward rounded rational bounds are the ones printed in the paper.
    paper_bounds={
      'F_low_interval':I(Q(48017,100000),Q(50674,100000)),
      'F_high_interval':I(Q(1005,10**10),Q(1165,10**10)),
      'DF_low_interval':I(-1683000,-1452000),
      'DF_high_interval':I(Q(-566,10**9),Q(-430,10**9))
    }
    for name,value in [('F_low_interval',image_low),('F_high_interval',image_high),
                       ('DF_low_interval',derivative_low),('DF_high_interval',derivative_high)]:
        bound=paper_bounds[name]
        require(bound.lo < value.lo <= value.hi < bound.hi)
    require(Q(1683000)*Q(566,10**9)==Q(952578,10**6)<Q(24,25))
    result={
      'status':'PASS',
      'proof':'Exact rational trapping intervals and Banach contraction certify an attracting two-cycle; a secondary IVT check certifies existence.',
      'parameters':{'alpha':'4/7','a21':'4/7','b1':'1/8','b2':'7/8','h':'1',
                    'A':[['-9990','10'],['9990','-10']],'stationary_first_component':'1/1000'},
      'left':str(left.lo),'right':str(right.lo),
      'F_left':left_step.record(),'F_right':right_step.record(),
      'F2_left':left_second.record(),'F2_right':right_second.record(),
      'F_bracket':whole_step.record(),
      'left_sign_margin_lower':str(left_second.lo-Q(1,4)),
      'right_sign_margin_lower':str(Q(3,4)-right_second.hi),
      'local_derivative_at_true_steady_state':str(local_R),
      'rational_initial_first_component':'1/2',
      'low_trapping_interval':low.record(),
      'high_trapping_interval':high.record(),
      'F_low_interval':image_low.record(),
      'F_high_interval':image_high.record(),
      'DF_low_interval':derivative_low.record(),
      'DF_high_interval':derivative_high.record(),
      'F2_derivative_upper_bound':str(contraction_bound),
      'F2_derivative_rational_bound':'24/25',
      'attracting_cycle_proof':'Banach contraction theorem on both closed trapping intervals, with disjoint images.',
      'paper_bounds':{name:bound.record() for name,bound in paper_bounds.items()},
      'certified_integer_fourth_root_bounds':ROOT_CHECKS,
      'floating_point_used_for_proof':False
    }
    out=Path(__file__).with_name('mprk22_cycle_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','proof','parameters','certified_integer_fourth_root_bounds','floating_point_used_for_proof']},indent=2))
    for name in ['F_left','F_right','F2_left','F2_right','F_bracket']:
        v=result[name]
        print(name, '[',v['decimal_lower_display_only'],',',v['decimal_upper_display_only'],']')
    for name in ['F_low_interval','F_high_interval','DF_low_interval','DF_high_interval']:
        v=result[name]
        print(name, '[',v['decimal_lower_display_only'],',',v['decimal_upper_display_only'],']')
    print('certified F^2 contraction bound < 24/25:',float(contraction_bound))
    print('local derivative:',local_R)

if __name__=='__main__': main()
