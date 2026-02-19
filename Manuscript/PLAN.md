# Blackout
## Research Outline
### Abstract
What explains conflict escalation? Past studies find different factors impacting conflict escalation ranging from structural factors to economic shocks, and organizational compositions. Recent technological innovations affecting escalation potential has often been neglected by these studies. By implementing new standards regarding conflict escalation, using DiD/DiDiD methods with noble dataset on the Arab Spring in Tunisia, we argue that internet blackout is likely to increase/decrease the likelihood of conflict escalation. This research contributes to past literature in three ways: 1) we integrate different measure of conflict escalation, which is expected to illustrate conflict escalation more specifically. 2) we show that fast changing technological environment can be a key aspect in conflict research nowadays, and 3) implement natural experiment to build causality.

### Research Question
What explains conflict escalation?

### Theory

Increase: Many people do not engage in the protest -> Internet blackout is universal for people -> Aggro to people -> More people, conflict type increase


Decrease:
1) Anticipation: 
    - Assumption: Internet blackout commonly precedes outright state repression
    - Internet blackout -> Signals potential participants -> Increased risks of getting repressed -> Deter participation -> Lowers probability of success in violent contention -> Less escalation
2) Disruption:
    - Assumption: Internet blackout cut out channel of mobilization through social media
    - Internet blackout -> Restricted mobilization -> Less participation -> Lowers probability of success in violent contention -> Less escalation

### Unit of Analysis
TBD (1st administrative unit-day/week or 2nd administrative unit-day/week)
- No year level; Maginot line is month

### Outcome
#### Definition
Conflict Escalation: Defined in three different ways
1) Tactical Escalation (Change in code of contention from nonviolence to violence)
    - Ordinal Coding
        - Protest to Protest
        - Protest to Riot
        - Protest to Terrorist attacks
        - Protest to Communal Violence (UCDP-GED 3) 
            -> Need to think about riot to communal violence
2) Numerical Escalation (Change in numbers that constitutes contentious political event)
    - Number of Participants
        - Ordinal Coding
    - Number of Fatalities/Casualties
        - Numerical Coding
        
        
#### Coding scheme
1) Find relevant events in ACLED
    - Retrack the original sources
    - Give the relevant events the same code - e.g., Tunisia_A_001
      - What are relevant events?
        1) Two or more events that are different in time demonstrate CLEAR common cause
        2) One or more events are different in time DIRECTLY state that they are inspired / influenced by precedent event(s)
        -> To YSL: Any further suggestions?
    - After finishing ACLED, go on to GTD and UCDP-GED '3'
    
2) Protest to Protest
    - Check the ordinal code and figure out whether the size increased or decreased (-1, 0, 1)
    - Check the casualty (if applicable) and fatality data
    
    
3) Protest to Riot / Terrorist attacks / Communal violence
    - If protest and riot, terrorist attacks, communal violence share the same code, then they are escalation! (No extra raw data checking)
    - Dichotomous coding (No - 0, Yes - 1)
    
#### Notes 
- to YSL: 
    - How can we deal with 'Spatial conflict escalation?' i.e., events in Adm-1 'a' influence conflict escalation in Adm-1 'b'?
        - A: 
            1) Can check spatial dependencies using different spatial models
                - Spatial Lag Model: Escalation in unit *it* influenced by escalation in unit *jt* (Spill over of outcome)
                - Spatial Lag of X Model: Treatment in unit *it* influenced by escalation in unit *jt* (Spillover of treatment)
                - Spatial Error Model: Estimated errors are spatially dependent between unit *it* and *jt*, rather than following normal distribution (i.e. being random)
                - Spatial Durbin Model: Generalized Spatial model that covers all models above
            2) Spatial models can be extended to panel data using FEs. There are computationally heavier options which are trickier than former using bayesian framework.
            3) In most cases, where the dependencies are weak, using TWFE themselves is acceptable.
    - Also, how can we deal with adm-1s w/o ANY protests?
        - A:
            - Two options to consider:
                1) Treat them as true zeroes:
                    - Code them as 0 (no escalation)
                    - We can preserve the number of observations and given escalation in nature is a rare event, I think it is theoretically more plausible
                    - Code them as NA and get rid of them from our observation risks: 1) we will lose out a lot of observations, 2) selection bias, and 3) distort spatial dependency structure.
                2) Potential for reporting bias can be addressed in robustness checks (Employable strategies need to be checked in later times)


### Treatment
Internet Penetration
-> Might benefit from binning (e.g, high-internet adm-1 / low-internet adm-1)
*
Internet Blackout
- To CHK: Wouldn't internet blackout be potentially endogenous to the level of internet penetration? If so, assumptions for DiD (i.e. parallel trends assumptions) might get violated. In that case, we will benefit from framing the study as estimating conditional effect rather than as DiD that strictly infer causality. Or is there any potential reasons that we can justify treating internet blackout exogeneous to escalation risk that varies locally by the level of internet penetration for our cases? 

### Covariate
GDPPC
Population
Distance from capital?


### Scope
Tunisia 2010 (Q4) - 2011 (Q1)

(Maybe Egypt in case Tunisian context is not applicable - e.g., VPN, proxy etc.)

### Method
TBD (DiD or DiDiD) w/ twfe

Tentative code: felm(escalation ~ Internet_penetration * Internet_shutdown + covariate | adm-1+year, data=data)

### Additional (Robustness) check
Observable implication
 - What is our alternative outcome instead of conflict escalation?
 - What is the condition at which our theory does not work?
 - What is the assumption are we making?
 
Placebo test (Lagging / adjusting internet blackout timing)

Pre-trend test (right before blackout, no difference between high and low internet adm-1)

To YSL: Additional Robustness check?



### Contribution
1) Newly defining the concept of escalation
2) Bridging digital repression literature with conflict escalation literature
3) Applying Natural experimental design for causal inference


## Timeline
### March 2026
#### Outcome
1) Check raw data of outcome (Tunisia 2010 Q4 - 2011 Q1)
    - ACLED & SCAD (event to campaign coding eligibility & source credibility)
2) Review dataset introduced by <strong><a href="https://journals.sagepub.com/doi/10.1177/0022002718777050" target="_blank" rel="noopener">Cunningham 2018</a></strong>

### April - June 2026
#### Outcome
1) Check raw data of outcome (Tunisia 2010 Q4 - 2011 Q1)
    - GTD & UCDP (event to campaign coding eligibility & source credibility)
2) Check raw data of outcome (Egypt)
    - ACLED & SCAD (event to campaign coding eligibility & source credibility)
#### Treatment
1) Blackout incident data collection
2) Internet penetration - NEED TO CHECK WHETHER THE DATA EXISTS
#### Lit Review
1) Escalation (YSL)
2) Blackout (CHK)
