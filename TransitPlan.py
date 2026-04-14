import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from astropy.visualization import astropy_mpl_style, quantity_support
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from IPython.display import HTML, display
from io import BytesIO
import base64
from dataclasses import dataclass

import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore', np.RankWarning)

from astroplan import Observer, FixedTarget ,Constraint, is_observable, observability_table
from astroplan import EclipsingSystem
from astroplan import (AltitudeConstraint, AirmassConstraint,
                       AtNightConstraint, MoonSeparationConstraint)
from astroplan.moon import moon_illumination
# from astropy.coordinates import get_moon
from astropy.coordinates import get_body
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive


def observing_night_from_sunset(observer, time):
    sunset = observer.sun_set_time(time, which="previous")
    local_dt = sunset.to_datetime(timezone=observer.timezone)
    return local_dt.strftime("%a. %Y-%m-%d")

def get_moon_data(observer, target, time):
    # Get moon illumination (0.0 - 1.0)
    # 0 = New Moon, 1 = Full Moon
    illumination = moon_illumination(time)
    
    # Get Moon coordinates from the observer's location
    moon_coord = get_body("moon", time, location=observer.location)
    
    # Calculate angular separation
    # separation = target.coord.separation(moon_coord)
    separation = moon_coord.separation(target.coord)
    
    return illumination, separation

def get_airmass(time,observer,target):
    return observer.altaz(time, target).secz

def alt_at(time,observer,target):
    return observer.altaz(time, target).alt

def box_transit(t, t_ing, t_egr, depth):
    y = np.zeros_like(t)
    in_transit = (t >= t_ing) & (t <= t_egr)
    y[in_transit] = -depth
    return y


def evaluate_transit(
    observer,
    target,
    ingress,
    mid,
    egress,
    constraints,
    baseline=None,
    nsamples=200
):
    duration = egress - ingress
    if baseline is None:
        baseline = 0.5 * duration
    else:
        baseline = baseline*u.hour

    t0 = ingress - baseline
    t1 = egress + baseline

    times = t0 + np.linspace(
        0, (t1 - t0).to_value(u.hour), nsamples
    ) * u.hour

    # Evaluate constraints
    observable = np.ones(len(times), dtype=bool)
    for c in constraints:
        observable &= c(observer, target, times)

    # Fractions
    in_transit = (times >= ingress) & (times <= egress)

    frac_transit = observable[in_transit].mean()
    frac_baseline = observable[~in_transit].mean()

    # Moon
    moon = get_body("moon",times, observer.location)
    moon_sep = moon.separation(target.coord).min()
    moon_illum = moon_illumination(mid)

    # Altitudes
    def alt(t):
        return observer.altaz(t, target).alt
    
    def airmass(t):
        return observer.altaz(t, target).secz
    
    airs= {"start": airmass(ingress),
           "mid": airmass(mid),
           "end": airmass(egress),
          }

    alts = {
        "start": alt(ingress),
        "mid": alt(mid),
        "end": alt(egress),
    }

    return {
        "times": times,
        "observable": observable,
        "frac_transit": frac_transit,
        "frac_baseline": frac_baseline,
        "moon_sep": moon_sep,
        "moon_illum": moon_illum,
        "alts": alts,
        "airs": airs,
    }

def transit_icon(times, ingress, egress, observable, depth):
    fig, ax = plt.subplots(figsize=(2.2, 0.6))

    y = np.zeros(len(times))
    in_transit = (times >= ingress) & (times <= egress)
    y[in_transit] = -depth

    # Non-observable (dashed)
    ax.plot(times.jd, y, color='0.7', lw=1, ls='--')

    # Observable (blue)
    ax.plot(times[observable].jd, y[observable], color='C0', lw=2)

    ax.set_axis_off()
    plt.tight_layout(pad=0.1)

    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=150)
    plt.close(fig)

    data = base64.b64encode(buf.getvalue()).decode()
    return f'<img src="data:image/png;base64,{data}"/>'

@dataclass
class TransitEvent:
    target: object
    ingress: Time
    mid: Time
    egress: Time
    depth: float
    mag: float


def create_transit_events(target, mid_transits, ingress_times, egress_times, params, constraints=None):
    """
    Create a list of TransitEvent objects from input parameters.
    
    Parameters:
    -----------
    target : FixedTarget
        The target object from astroplan
    mid_transits : list or array of Time objects
        Mid-transit times
    ingress_times : list or array of Time objects
        Ingress times
    egress_times : list or array of Time objects
        Egress times
    params : dict
        Dictionary containing 'pl_trandep' (transit depth) and 'sy_vmag' (V magnitude)
    constraints : list, optional
        List of observability constraints
        
    Returns:
    --------
    events : list
        List of TransitEvent objects
    """
    events = []
    
    for mid, ing, egr in zip(mid_transits, ingress_times, egress_times):
        events.append(
            TransitEvent(
                target=target,
                ingress=ing,
                mid=mid,
                egress=egr,
                depth=np.array(params['pl_trandep']),
                mag=np.array(params['sy_vmag']),
            )
        )
    
    return events

def display_transits(target, mid_transits, ingress_times, egress_times, params, observer, constraints,baseline=None, visibility_cut=None,plot=True):
    """
    Generate a DataFrame summarizing transit observability for multiple events given target and observatory objects from astroplan.

    Parameters:
    -----------
    target : FixedTarget
        The target object from astroplan
    mid_transits : list or array of Time objects
        Mid-transit times
    ingress_times : list or array of Time objects
        Ingress times
    egress_times : list or array of Time objects
        Egress times
    params : dict
        Exoplanet archive dictionary containing 'pl_trandep' (transit depth) and 'sy_vmag' (V magnitude)
    observer : Observeratory object
        The observatory object from astroplan
    constraints : list
        List of observability constraints from astroplan
    baseline : float (time), optional
        Baseline duration before and after transit for observability calculations 
        Should be baseline time in hours desired before and after transit, not the total available baseline
        defaults to 0.5*transit duration if None
    visibility_cut : float, optional
        Percentage of transit visibility to include in output DataFrame
    plot : bool, optional
        If True, display the DataFrame as an HTML table. If False, return the DataFrame.
    Returns:
    --------
    df : pandas DataFrame
        DataFrame summarizing transit observability metrics for each event

    """
    events = create_transit_events(target, mid_transits, ingress_times, egress_times, params, constraints)
    
    rows = []
    
    for event in events:
        metrics = evaluate_transit(
            observer,
            event.target,
            event.ingress,
            event.mid,
            event.egress,
            constraints,
            baseline=baseline,
            nsamples=200
        )
        
        icon = transit_icon(
            metrics["times"],
            event.ingress,
            event.egress,
            metrics["observable"],
            event.depth
        )
        
        rows.append({
            "Local Evening Date": f"{observer.name}<br>"
            f"{observing_night_from_sunset(observer, event.mid)}<br>"
            f"({observer.timezone})",
            "Moon": (
                f"{metrics['moon_illum']:.0%} @ "
                f"{metrics['moon_sep'].to(u.deg).value:.0f}°"
            ),
            "Start-Mid-End": (
                f"{event.ingress.strftime('%H:%M')}<br>"
                f"{event.mid.strftime('%H:%M')}<br>"
                f"{event.egress.strftime('%H:%M')}"
            ),
            "Alt (S/M/E)": (
                f"{metrics['alts']['start'].value:.0f}°<br>"
                f"{metrics['alts']['mid'].value:.0f}°<br>"
                f"{metrics['alts']['end'].value:.0f}°"
            ),
            "Airmass (S/M/E)": (
                f"{metrics['airs']['start'].value:.2f}<br>"
                f"{metrics['airs']['mid'].value:.2f}<br>"
                f"{metrics['airs']['end'].value:.2f}"
            ),
            "Transit": icon,
            "% Transit": 100*metrics['frac_transit'],
            "% Baseline": 100*metrics['frac_baseline'],
            "V mag": float(event.mag),
            "Depth (ppt)": float(event.depth)*10,
            "Target": event.target.name,
        })
    
    
    df = pd.DataFrame(rows)
    if visibility_cut is not None:
        df = df[df["% Transit"] >= visibility_cut]

    if plot:
        return (df.style
                .format({'% Transit':'{:.0f}%','% Baseline':'{:.0f}%','Depth (ppt)':'{:.2f}','V mag':'{:.2f}'})
                .set_properties(**{'background-color':"#FFFFFF",'color':'black'})
                .set_table_styles([{"selector":"th.col_heading","props":"text-align:center; font-size:1.5em; background-color: #000066;"}])
                .background_gradient(subset=['% Transit','% Baseline'],cmap="RdYlGn",vmin=0,vmax=100))
        # return display(HTML(df.style.set_properties(**{'background-color': "#FFFFFF", 'color': 'black'}).to_html(escape=False, index=False)))
    else:
        return df


def simultaneous_observing(target, mid_transits, ingress_times, egress_times, params, observatory_1,comparison_observatory, constraints, baseline=None, visibility_cut=None,plot=True):
    """
    Compare transit observability across multiple observatories for a given target and transit events.

    Parameters:
    -----------
    target : FixedTarget
        The target object from astroplan
    mid_transits : list or array of Time objects
        Mid-transit times
    ingress_times : list or array of Time objects
        Ingress times
    egress_times : list or array of Time objects
        Egress times
    params : dict
        Exoplanet archive dictionary containing 'pl_trandep' (transit depth) and 'sy_vmag' (V magnitude)
    observatory_1 : Observatory object from astroplan
        First observatory for comparison
    comparison_observatory : Observatory object from astroplan
        Second observatory for comparison
    constraints : list
        List of observability constraints from astroplan
    baseline : float (time), optional
        Baseline duration before and after transit for observability calculations 
        Should be baseline time in hours desired before and after transit, not the total available baseline
        defaults to 0.5*transit duration if None
    visibility_cut : float, optional
        Percentage of transit visibility to include in output DataFrame
    plot : bool, optional
        Whether to return a html styled DataFrame for plotting (default is True)
    Returns:
    --------
    df : pandas DataFrame
        DataFrame summarizing transit observability metrics for each event and observatory

    """
    all_rows = []
    
    
    events = create_transit_events(target, mid_transits, ingress_times, egress_times, params, constraints)
    
    for event in events:
        metrics = evaluate_transit(
            observatory_1,
            event.target,
            event.ingress,
            event.mid,
            event.egress,
            constraints,
            baseline=baseline,
            nsamples=200
        )
        icon = transit_icon(
            metrics["times"],
            event.ingress,
            event.egress,
            metrics["observable"],
            event.depth
        )

        metrics2 = evaluate_transit(
            comparison_observatory,
            event.target,
            event.ingress,
            event.mid,
            event.egress,
            constraints,
            baseline=baseline,
            nsamples=200
        )
        icon2 = transit_icon(
            metrics2["times"],
            event.ingress,
            event.egress,
            metrics2["observable"],
            event.depth
        )

        all_rows.append({
            f"Local Evening Date @<br> {observatory_1.name}": f"{observatory_1.name}<br>"
            f"{observing_night_from_sunset(observatory_1, event.mid)}<br>"
            f"({observatory_1.timezone})",
            "Moon": (
                f"{metrics['moon_illum']:.0%} @ "
                f"{metrics['moon_sep'].to(u.deg).value:.0f}°"
            ),
            "Start-Mid-End": (
                f"{event.ingress.strftime('%H:%M')}<br>"
                f"{event.mid.strftime('%H:%M')}<br>"
                f"{event.egress.strftime('%H:%M')}"
            ),
            "Alt (S/M/E)": (
                f"{metrics['alts']['start'].value:.0f}°<br>"
                f"{metrics['alts']['mid'].value:.0f}°<br>"
                f"{metrics['alts']['end'].value:.0f}°"
            ),
            "Airmass (S/M/E)": (
                f"{metrics['airs']['start'].value:.2f}<br>"
                f"{metrics['airs']['mid'].value:.2f}<br>"
                f"{metrics['airs']['end'].value:.2f}"
            ),
            f"{event.target.name} Transit": icon,
            
            "% Transit": 100*metrics['frac_transit'],
            "% Baseline": 100*metrics['frac_baseline'],
            f"Local Evening Date @<br> {comparison_observatory.name}": f"{comparison_observatory.name}<br>"
            f"{observing_night_from_sunset(comparison_observatory, event.mid)}<br>"
            f"({comparison_observatory.timezone})",
            f"Transit @ {comparison_observatory.name}": icon2,
            f"Moon @ {comparison_observatory.name}": (
                f"{metrics2['moon_illum']:.0%} @ "
                f"{metrics2['moon_sep'].to(u.deg).value:.0f}°"
            ),
            # "% Transit": 100*metrics['frac_transit'],
            # "% Baseline": 100*metrics['frac_baseline'],
            f"% Transit @<br> {comparison_observatory.name}": 100*metrics2['frac_transit'],
            f"% Baseline @<br> {comparison_observatory.name}": 100*metrics2['frac_baseline'],
        })

    df = pd.DataFrame(all_rows)
    if visibility_cut is not None:
        df = df[df["% Transit"] >= visibility_cut]
    if plot:
        return (df.style
                .format({'% Transit':'{:.0f}%','% Baseline':'{:.0f}%','Depth (ppt)':'{:.2f}',f"% Transit @<br> {comparison_observatory.name}":'{:.0f}%',f"% Baseline @<br> {comparison_observatory.name}":'{:.0f}%',f"% Transit @<br> {comparison_observatory.name}":'{:.0f}%'})
                .set_properties(**{'background-color':"#FFFFFF",'color':'black',})
                .set_table_styles([{"selector":"th.col_heading","props":"text-align:center; font-size:1.5em; background-color: #000066;"}])
                .background_gradient(subset=['% Transit','% Baseline'],cmap="RdYlGn",vmin=0,vmax=100)
                .background_gradient(subset=[f"% Transit @<br> {comparison_observatory.name}",f"% Baseline @<br> {comparison_observatory.name}"],cmap="RdYlGn",vmin=0,vmax=100))
    else:
        return df

def single_transit(target, mid_transits, ingress_times, egress_times, params, observer, constraints, baseline=None, visibility_cut=None,plot=True,event_index=0):
    """
    Generate a DataFrame summarizing transit observability for single event determined using either display_transits or simultaneous_observing functions, given target and observatory objects from astroplan.

    Parameters:
    -----------
    target : FixedTarget
        The target object from astroplan
    mid_transits : list or array of Time objects
        Mid-transit times
    ingress_times : list or array of Time objects
        Ingress times
    egress_times : list or array of Time objects
        Egress times
    params : dict
        Exoplanet archive dictionary containing 'pl_trandep' (transit depth) and 'sy_vmag' (V magnitude)
    observer : Observeratory object
        The observatory object from astroplan
    constraints : list
        List of observability constraints from astroplan
    baseline : float (time), optional
        Baseline duration before and after transit for observability calculations 
        Should be baseline time in hours desired before and after transit, not the total available baseline
        defaults to 0.5*transit duration if None
    visibility_cut : float, optional
        Percentage of transit visibility to include in output DataFrame
    plot : bool, optional
        If True, display the DataFrame as an HTML table. If False, return the DataFrame.
    event_index : int
        Index of the transit event to display (default is 0, but depending on constraints and visibility_cut imposed, this will likely not be 0)
        Determine using the output of display_transits or simultaneous_observing 
    Returns:
    --------
    df : pandas DataFrame
        DataFrame summarizing transit observability metrics for single event with additional metrics at the baseline start and end times

    """
    events = create_transit_events(target, mid_transits, ingress_times, egress_times, params, constraints)
    print(events)
    event = events[event_index]

    rows = []

    # for event in events:
    metrics = evaluate_transit(
        observer,
        event.target,
        event.ingress,
        event.mid,
        event.egress,
        constraints,
        baseline=baseline,
        nsamples=200
    )
    
    icon = transit_icon(
        metrics["times"],
        event.ingress,
        event.egress,
        metrics["observable"],
        event.depth
    )

    baseline_start = event.ingress - baseline*u.hour
    baseline_end = event.egress + baseline*u.hour
    bair_start = get_airmass(baseline_start, observer, event.target)
    bair_end = get_airmass(baseline_end, observer, event.target)
    balt_start = alt_at(baseline_start, observer, event.target)
    balt_end = alt_at(baseline_end, observer, event.target)
    
    rows.append({
            "Local Evening Date": f"{observer.name}<br>"
            f"{observing_night_from_sunset(observer, event.mid)}<br>"
            f"({observer.timezone})",
            "Moon": (
                f"{metrics['moon_illum']:.0%} @ "
                f"{metrics['moon_sep'].to(u.deg).value:.0f}°"
            ),
            "Start-Mid-End": (
                f"<span style='opacity:0.6'>{baseline_start.strftime('%H:%M')}</span><br>"
                f"{event.ingress.strftime('%H:%M')}<br>"
                f"{event.mid.strftime('%H:%M')}<br>"
                f"{event.egress.strftime('%H:%M')}<br>"
                f"<span style='opacity:0.6'>{baseline_end.strftime('%H:%M')}</span>"
            ),
            "Alt (S/M/E)": (
                f"<span style='opacity:0.6'>{balt_start.value:.0f}°</span><br>"
                f"{metrics['alts']['start'].value:.0f}°<br>"
                f"{metrics['alts']['mid'].value:.0f}°<br>"
                f"{metrics['alts']['end'].value:.0f}°<br>"
                f"<span style='opacity:0.6'>{balt_end.value:.0f}°</span>"
            ),
            "Airmass (S/M/E)": (
                f"<span style='opacity:0.6'>{bair_start.value:.2f}</span><br>"
                f"{metrics['airs']['start'].value:.2f}<br>"
                f"{metrics['airs']['mid'].value:.2f}<br>"
                f"{metrics['airs']['end'].value:.2f}<br>"
                f"<span style='opacity:0.6'>{bair_end.value:.2f}</span>"
            ),
            "Transit": icon,
            "% Transit": 100*metrics['frac_transit'],
            "% Baseline": 100*metrics['frac_baseline'],
            "V mag": float(event.mag),
            "Depth (ppt)": float(event.depth)*10,
            "Target": event.target.name,
        })


    df = pd.DataFrame(rows)
    if visibility_cut is not None:
        df = df[df["% Transit"] >= visibility_cut]

    if plot:
        return (df.style
                .format({'% Transit':'{:.0f}%','% Baseline':'{:.0f}%','Depth (ppt)':'{:.2f}','V mag':'{:.2f}'})
                .set_properties(**{'background-color':"#FFFFFF",'color':'black'})
                .set_table_styles([{"selector":"th.col_heading","props":"text-align:center; font-size:1.5em; background-color: #000066;"}])
                .background_gradient(subset=['% Transit','% Baseline'],cmap="RdYlGn",vmin=0,vmax=100))
        # return display(HTML(df.style.set_properties(**{'background-color': "#FFFFFF", 'color': 'black'}).to_html(escape=False, index=False)))
    else:
        return df
    # return display_transits(target, mid_transits, ingress_times, egress_times, params, observer, constraints, baseline=baseline, visibility_cut=visibility_cut)