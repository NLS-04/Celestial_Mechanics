using System;
using static System.Math;
using Vectors;
using System.Collections.Generic;

namespace Celestial_Mechanics {
	public class Constants {
		public const double G   = 6.67430E-11d;
		public const double e   = Math.E;
		public const double pi  = Math.PI;
		public const double tau = 2 * pi;

		/// <summary>convert a degrees value into a radians value by multiplying with this constant <see cref="rad"/></summary>
		public const double rad = pi / 180;

		/// <summary>convert a radians value into a degrees value by multiplying with this constant <see cref="deg"/></summary>
		public const double deg = 180 / pi;
	}

	public class Body {
		// mu = G*M
		/// <summary>Standard gravitational parameter 'mu' of an <see cref="Body"/> in [m^3/s^2]</summary>
		public readonly double mu;
		/// <summary>Mass of an <see cref="Body"/> in [kg]</summary>
		public readonly double mass;
		/// <summary>Average radial distance of the surface of an <see cref="Body"/> in [m]</summary>
		public readonly double radius;

		public readonly Direction orientation = Direction.defaultDirection;


		/// <summary>Constructor for an (large) orbitable Body</summary>
		/// <param name="mu">Standard gravitational parameter 'mu' of <see langword="this"/> <see cref="Body"/> in [m^3/s^2]</param>
		/// <param name="radius">Average radial distance of the surface of <see langword="this"/> <see cref="Body"/> in [m]</param>
		/// <param name="orient">A struct to orientate a coordinate system of <see langword="this"/> <see cref="Body"/>; if its <see langword="null"/> the <see cref="Body.orientation"/> is set to the <see cref="defaultDirection"/></param>
		public Body( float mu, double radius /*, Direction? orient=null*/ ) {
			this.mu     = mu;
			this.mass   = mu / Constants.G;
			this.radius = radius;
			//this.orientation = orient ?? Direction.defaultDirection;
		}

		/// <summary>Constructor for an (large) orbitable Body</summary>
		/// <param name="radius">Average radial distance of the surface of <see langword="this"/> <see cref="Body"/> in [m]</param>
		/// <param name="mass">Mass of <see langword="this"/> <see cref="Body"/> in [kg]</param>
		/// <param name="orient">A struct to orientate a coordinate system of <see langword="this"/> <see cref="Body"/>; if its <see langword="null"/> the <see cref="Body.orientation"/> is set to the <see cref="defaultDirection"/></param>
		public Body( double radius, float mass /*, Direction? orient=null*/ ) {
			this.mass   = mass;
			this.mu     = mass * Constants.G;
			this.radius = radius;
			//this.orientation = orient ?? Direction.defaultDirection;
		}

		public readonly static Body KERBIN = new Body( 3.5316E12f,      600_000d   );
		public readonly static Body EARTH  = new Body( 3.986004418E14f, 6_371_000d );
		public readonly static Body MARS   = new Body( 4.282837E13f,	3_389_500d );

		public override bool Equals( object obj ) {
			if ( !(obj is Body) )
				return false;
			Body bod = (Body) obj;
			return mu == bod.mu && radius == bod.radius && ( orientation.Equals(bod.orientation) );
		}

		public static bool operator ==( Body a, Body b ) => a.Equals( b );
		public static bool operator !=( Body a, Body b ) =>!a.Equals( b );
	}

	/// <summary>
	/// 
	/// </summary>
	public class Orbit {
		/* LAN (longitude of Ascending Node) measured to X-Axis of the Body
		 * inclination measured to XY-Plane of the Body
		 * Body's North == Z-Axis
		 */

		/* All (readable) angle Values are in [DEGREES] but are internaly converted and calculated as [RADIANS] */

		public Body body;

		// vectors are normalized, values are stored at their associativ double variable
		public readonly Vector eccentricity_Vector;
		public readonly Vector angularMomentum_Vector;
		public readonly Vector nodeLine_Vector;

		public readonly double apoapsis, periapsis, semiMajorAxis, semiMinorAxis;
		public readonly double inclination, eccentricity, angularMomentum, longitudeOfAscendingNode, argumentOfPeriapsis;
		public readonly double period, epoch, meanAnomaly_At_Epoch, semi_latus_rectum;

		public ManeuverNode? maneuverNode { get; set; }

		private readonly Direction dir_at_periapsis;

		#region constructors
		/// <summary><see cref="Orbit"/> constructor for general Orbits, orientated in 3D</summary>
		/// <param name="inc">Inclination in [degrees]</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <param name="sma">Semi-Major-Axis of Orbit in [m]</param>
		/// <param name="lan">Longitude of Ascending Node in [degrees]</param>
		/// <param name="argOfPer">Argument of Periapsis in [degrees]</param>
		/// <param name="tEpoch">Arbitrary timestamp to calculate time dependend values of this orbit (e.g. position, velocity etc.)</param>
		/// <param name="body"><see cref="Body"/> to be orbited around</param>
		public Orbit( double inc, double ecc, double sma, double lan, double argOfPer, double mEp, double tEpoch, Body body ) {
			inclination              = inc;
			eccentricity             = ecc;
			semiMajorAxis            = sma;
			longitudeOfAscendingNode = lan;
			argumentOfPeriapsis      = argOfPer;
			epoch                    = tEpoch;
			this.body				 = body;

			apoapsis  = sma * ( 1 + ecc );
			periapsis = sma * ( 1 - ecc );

			semiMinorAxis   = calculate_semiMinorAxis( semiMajorAxis, eccentricity );
			angularMomentum = calculate_angularMomentum( semiMajorAxis, eccentricity, body.mu );
			period          = calculate_Period( semiMajorAxis, body.mu );

			semi_latus_rectum = angularMomentum * angularMomentum / body.mu;


			nodeLine_Vector        = Vector.angleAxis( body.orientation.forward, body.orientation.top, longitudeOfAscendingNode * Constants.rad ).normalized();
			angularMomentum_Vector = Vector.angleAxis( body.orientation.top, nodeLine_Vector, inclination * Constants.rad ).normalized();
			eccentricity_Vector    = Vector.angleAxis( nodeLine_Vector, angularMomentum_Vector, argOfPer * Constants.rad ).normalized();

			meanAnomaly_At_Epoch = mEp;

			dir_at_periapsis = Direction.lookDirUp( eccentricity_Vector, angularMomentum_Vector );
		}

		/// <summary><see cref="Orbit"/> constructor for general Orbits, orientated in 3D, based on an initial position and velocity</summary>
		/// <param name="pos">initial position of satellite in [m]</param>
		/// <param name="vel">initial velocity of satellite in [m/s]</param>
		/// <param name="tEpoch">Arbitrary timestamp to calculate time dependend parameters of this orbit (e.g. position, velocity, flight angle etc.)</param>
		/// <param name="body"><see cref="Body"/> to be orbited around</param>
		public Orbit( Vector pos, Vector vel, double tEpoch, Body body ) {
			epoch     = tEpoch;
			this.body = body;

			angularMomentum_Vector = Vector.cross( pos, vel );
			angularMomentum = angularMomentum_Vector.mag();

			inclination = Vector.vang(angularMomentum_Vector, body.orientation.top) * Constants.deg;

			nodeLine_Vector = Vector.cross( body.orientation.top, angularMomentum_Vector );

			longitudeOfAscendingNode = Acos( nodeLine_Vector.x / nodeLine_Vector.mag() ) * Constants.deg;
			longitudeOfAscendingNode = nodeLine_Vector.y < 0 ? 360 - longitudeOfAscendingNode : longitudeOfAscendingNode;

			eccentricity_Vector = ( 1 / body.mu ) * ( Vector.cross( vel, angularMomentum_Vector ) - body.mu * pos.normalized() );
			eccentricity = eccentricity_Vector.mag();

			argumentOfPeriapsis = Vector.vang( nodeLine_Vector, eccentricity_Vector ) * Constants.deg;
			argumentOfPeriapsis = eccentricity_Vector.z < 0 ? 360 - argumentOfPeriapsis : argumentOfPeriapsis;

			apoapsis  = altitude( 180 );
			periapsis = altitude(   0 );

			semiMajorAxis = .5f * ( apoapsis + periapsis );
			semiMinorAxis = calculate_semiMinorAxis( semiMajorAxis, eccentricity );
			period = calculate_Period( semiMajorAxis, body.mu );

			semi_latus_rectum = angularMomentum * angularMomentum / body.mu;

			nodeLine_Vector = nodeLine_Vector.normalized();
			angularMomentum_Vector = angularMomentum_Vector.normalized();
			eccentricity_Vector = eccentricity_Vector.normalized();

			meanAnomaly_At_Epoch = meanAnomaly( Vector.vang( body.orientation.forward, eccentricity_Vector ) * Constants.deg );

			dir_at_periapsis = Direction.lookDirUp( eccentricity_Vector, angularMomentum_Vector );
		}
		#endregion

		#region operators
		public static bool operator ==( Orbit a, Orbit b ) => a.Equals( b );
		public static bool operator !=( Orbit a, Orbit b ) =>!a.Equals( b );

		public override bool Equals( object obj ) {
			if ( !(obj is Orbit) )
				return false;
			Orbit b = (Orbit) obj;
			return
				this.inclination == b.inclination &&
				this.eccentricity == b.eccentricity &&
				this.semiMajorAxis == b.semiMajorAxis &&
				this.longitudeOfAscendingNode == b.longitudeOfAscendingNode &&
				this.argumentOfPeriapsis == b.argumentOfPeriapsis &&
				this.body == b.body;
		}
		#endregion

		#region static_methods
		/// <summary>Calculate the time an satellite takes to make one revolution on an <see cref="Orbit"/></summary>
		/// <param name="sma">Semi-Major-Axis of Orbit in [meters]</param>
		/// <param name="mu">Standard gravitational parameter 'mu'(<seealso cref="Body.mu"/>) of an <see cref="Body"/> in [meters^3/seconds^2]</param>
		/// <returns>Orbits Period in [seconds]</returns>
		public static double calculate_Period( double sma, double mu ) {
			return 2 * Constants.pi * Sqrt( Pow( sma, 3 ) / mu );
		}

		/// <summary>Calculate the Angular Momentum</summary>
		/// <param name="sma">Semi-Major-Axis of Orbit in [meters]</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <param name="mu">Standard gravitational parameter 'mu'(<seealso cref="Body.mu"/>) of an <see cref="Body"/> in [meters^3/seconds^2]</param>
		/// <returns>Angular momentum in [meters^2/seconds]</returns>
		public static double calculate_angularMomentum( double sma, double ecc, double mu ) {
			return Sqrt( sma * mu * ( 1 - ecc * ecc ) );
		}

		/// <summary>Calculate the Semi Minor Axis of an <see cref="Orbit"/></summary>
		/// <param name="sma">Semi-Major-Axis of Orbit in [meters]</param>
		/// <returns>Semi-Minor-Axis in [meters]</returns>
		public static double calculate_semiMinorAxis( double sma, double ecc ) {
			return sma * Sqrt( 1 - ecc * ecc );
		}


		/// <summary>Calculate the radial distance of a satellite at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <param name="h">Angular momentum in [meters^2/seconds]</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <param name="mu">Standard gravitational parameter 'mu'(<seealso cref="Body.mu"/>) of an <see cref="Body"/> in [meters^3/seconds^2]</param>
		/// <returns>The radial distance in [seconds]</returns>
		public static double altitude( double trueAnomaly, double h, double ecc, double mu ) {
			return h * h / ( mu * ( 1 + ecc * Cos( trueAnomaly * Constants.rad ) ) );
		}

		/// <summary>Calculate the angle of the velocity Vector to the local horzion Plane at a True Anomaly</summary>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>Flight angle in [Degrees]</returns>
		public static double flight_angle( double ecc, double trueAnomaly ) {
			trueAnomaly *= Constants.rad;
			return Atan( ecc * Sin( trueAnomaly ) / ( 1 + ecc * Cos( trueAnomaly ) ) ) * Constants.deg;
		}

		#region Orbit position determination
		/// <summary>Calculate the Eccentricity Anomaly for a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <returns>The Eccentricity Anomaly in !! [radians] !!</returns>
		public static double eccentricityAnomaly( double trueAnomaly, double ecc ) {
			// the Eccentricity Anomaly in [radians]
			return 2.0d * Atan( Sqrt( ( 1.0d - ecc ) / ( 1.0d + ecc ) ) * Tan( .5d * trueAnomaly * Constants.rad ) );
		}

		/// <summary>Calculate the Mean Anomaly for a given Eccentricity Anomaly</summary>
		/// <param name="E">The Eccentricity Anomaly in [radians]</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <returns>The Mean Anomaly in [degrees]</returns>
		public static double meanAnomaly_from_eccentricityAnomaly( double E, double ecc ) {
			return (E - ecc * Sin( E )) * Constants.deg; // E in radians, returns in deg
		}

		/// <summary>Calculate the Mean Anomaly for a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <returns>The Mean Anomaly in [degrees]</returns>
		public static double meanAnomaly( double trueAnomaly, double ecc ) {
			return meanAnomaly_from_eccentricityAnomaly( eccentricityAnomaly( trueAnomaly, ecc ), ecc );
		}


		/// <summary>Calculate the time since periapsis for a given Mean Anomaly</summary>
		/// <param name="meanAnomaly">The Mean Anomaly in [degrees]</param>
		/// <param name="period">The orbits Period time in [seconds]</param>
		/// <returns>The time that passed since periapsis passage in [seconds]</returns>
		public static double meanAnomaly_to_time( double meanAnomaly, double period ) {
			return period * meanAnomaly*Constants.rad / Constants.tau;
		}

		/// <summary>Calculate the time since periapsis for a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <param name="period">The orbits Period time in [seconds]</param>
		/// <returns>The time that passed since periapsis passage in [seconds]</returns>
		public static double trueAnomaly_to_time( double trueAnomaly, double ecc, double period ) {
			return meanAnomaly_to_time( meanAnomaly( trueAnomaly, ecc ), period );
		}
		#endregion
		
		#endregion

		#region methods
		/// <summary>Calculate the radial distance of a satellite at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The radial distance for the given True Anomaly</returns>
		public double altitude( double trueAnomaly ) {
			return altitude( trueAnomaly, angularMomentum, eccentricity, body.mu );
		}

		/// <summary>Calculate the angle of the velocity to the local horzion at a True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>Flight angle in [Degrees]</returns>
		public double flight_angle( double trueAnomaly ) {
			return Orbit.flight_angle( this.eccentricity, trueAnomaly );
		}


		/// <summary>Calculate the velocity tangential to the Flight Path</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>Tangential velocity in [meters/second]</returns>
		public double tangential_velocity( double trueAnomaly ) {
			return angularMomentum / altitude( trueAnomaly );
		}

		/// <summary>Calculate the velocity radial outwards to the Flight Path</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>Radial outwards velocity in [meters/second]</returns>
		public double radial_velocity( double trueAnomaly ) {
			return body.mu / angularMomentum * eccentricity * Sin(trueAnomaly);
		}

		/// <summary>Calculate the velocity Tuple (tangetial, radial outwards) to the Flight Path</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>(tangetial velocity, radial outwards velocity) in [meters/second]</returns>
		public (double, double) velocities_tangential_radial( double trueAnomaly ) {
			return (tangential_velocity( trueAnomaly ), radial_velocity( trueAnomaly ));
		}


		#region Orbit position determination
		/// <summary>Calculate the Eccentricity Anomaly for a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The Eccentricity Anomaly in !!! [radians] !!</returns>
		public double eccentricityAnomaly( double trueAnomaly ) {
			return eccentricityAnomaly( trueAnomaly, eccentricity );
		}

		/// <summary>Calculate the Mean Anomaly for a given Eccentricity Anomaly</summary>
		/// <param name="E">The Eccentricity Anomaly in [radians]</param>
		/// <returns>The Mean Anomaly in [degrees]</returns>
		public double meanAnomaly_from_eccentricityAnomaly( double E ) {
			return meanAnomaly_from_eccentricityAnomaly( E, eccentricity );
		}

		/// <summary>Calculate the Mean Anomaly for a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The Mean Anomaly in [degrees]</returns>
		public double meanAnomaly( double trueAnomaly ) {
			return meanAnomaly( trueAnomaly, eccentricity );
		}


		/// <summary>Calculate the time since periapsis for a given Mean Anomaly</summary>
		/// <param name="meanAnomaly">The Mean Anomaly in [degrees]</param>
		/// <returns>The time that passed since periapsis passage in [seconds]</returns>
		public double meanAnomaly_to_time( double meanAnomaly ) {
			return meanAnomaly_to_time( meanAnomaly, period );
		}

		/// <summary>Calculate the time since periapsis for a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The time that passed since periapsis passage in [seconds]</returns>
		public double trueAnomaly_to_time( double trueAnomaly ) {
			return trueAnomaly_to_time( trueAnomaly, eccentricity, period );
		}

		/// <summary>Calculate the Mean Anomaly at a given time</summary>
		/// <param name="time">Time in [seconds]</param>
		/// <param name="TM">The Time mode controls how the time value should be interpreted</param>
		/// <returns>The Mean Anomaly in [degrees]</returns>
		public double meanAnomaly_at_time( double time, TimeMode TM=TimeMode.UT ) {
			double Me = 0d;
			if (TM == TimeMode.UT )
				Me = Constants.tau * ( time - epoch ) / period + meanAnomaly_At_Epoch;
			else 
				Me = Constants.tau * time / period;

			return ( Me % Constants.tau ) * Constants.deg;
		}

		/// <summary>Calculate the True Anomaly at a given time; the solution is non trivial, thus this value gets calculated via an iterative Algorithm</summary>
		/// <param name="time">Time in [seconds]</param>
		/// <param name="TM">The Time mode controls how the time value should be interpreted</param>
		/// <returns>The True Anomaly in [degrees]</returns>
		public double trueAnomaly_at_time( double time, TimeMode TM=TimeMode.UT ) {
			double Me = meanAnomaly_at_time( time, TM ) * Constants.rad;

			Solver.ERROR_THRESHOLD = 1E-5f;
			Solver.MAX_ITERATION_STEPS = 1_000u;

			Delegate func  = Solver.get_EccAnomaly_method( Solver.EccAnomalyMode.newton_LIT_n5 );
			var methodinfo = func.Method;

			return (double) methodinfo.Invoke( null, new object[] { Me, eccentricity } ) * Constants.deg;
		}
		#endregion


		/// <summary>Calculate the Direction at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		public Direction direction_at_trueAnomaly( double trueAnomaly ) {
			return dir_at_periapsis.angleAxis( angularMomentum_Vector, trueAnomaly );
		}

		#region State Vector methods
		/// <summary>Calculate the normalized Position Vector (e.g. the direction) at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The Position 3d-normalized-Vector </returns>
		public Vector normalized_position( double trueAnomaly ) {
			// r_hat or u_hat, position direction
			return Vector.angleAxis( eccentricity_Vector, angularMomentum_Vector, trueAnomaly ).normalized();
		}

		/// <summary>Calculate the Position Vector at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The Position 3d-Vector in [meters]</returns>
		public Vector position_at_trueAnomaly( double trueAnomaly ) {
			return altitude( trueAnomaly ) * normalized_position( trueAnomaly );
		}

		/// <summary>Calculate the Velocity Vector at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The Velocity 3d-Vector [meters/second]</returns>
		public Vector velocity_at_trueAnomaly( double trueAnomaly ) {
			Direction dir = direction_at_trueAnomaly( trueAnomaly );
			Vector radial_out_vel = radial_velocity( trueAnomaly ) * dir.top;
			Vector tangential_vel = tangential_velocity( trueAnomaly ) * Vector.cross( angularMomentum_Vector, dir.forward );
			return radial_out_vel + tangential_vel;
		}

		/// <summary>Calculate the Position Vector at a given Time</summary>
		/// <param name="time">Time in [seconds]</param>
		/// <param name="TM">The Time mode controls how the time value should be interpreted</param>
		/// <returns>The Position 3d-Vector in [meters]</returns>
		public Vector position_at_time( double time, TimeMode TM=TimeMode.UT ) {
			return position_at_trueAnomaly( trueAnomaly_at_time( time, TM ) );
		}

		/// <summary>Calculate the Velocity Vector at a given Time</summary>
		/// <param name="time">Time in [seconds]</param>
		/// <param name="TM">The Time mode controls how the time value should be interpreted</param>
		/// <returns>The Velocity 3d-Vector [meters/second]</returns>
		public Vector velocity_at_time( double time, TimeMode TM = TimeMode.UT ) {
			return velocity_at_trueAnomaly( trueAnomaly_at_time( time, TM ) );
		}
		#endregion

		public string staticInformation() {
			return $@"apoapsis:	{apoapsis}
periapsis:	{periapsis}
SMA:		{semiMajorAxis}
SMinorA:	{semiMinorAxis}
inc:		{inclination}
ecc:		{eccentricity}
h:		{angularMomentum}
LAN:		{longitudeOfAscendingNode}
argPer:		{argumentOfPeriapsis}
period:		{period}

epoch:		{epoch}
meanEp:		{meanAnomaly_At_Epoch}

ecc_Vec:	{eccentricity_Vector}
h_Vec:		{angularMomentum_Vector}
node_Vec:	{nodeLine_Vector}";
		}
		#endregion

		public enum TimeMode {
			/// <summary>Universal Time is the time in [seconds] that passed since an arbitrary time stamp</summary>
			UT,
			/// <summary>PeriapsRelativ is the time in [seconds] that passed since the Periapsis and is interpreted as 0 at Periapsis</summary>
			PeriapsRelativ
		}
	}

	public struct Direction {
		public readonly static Direction defaultDirection = new Direction( Vector.forward, Vector.side, Vector.up );

		public readonly Vector forward;
		public readonly Vector top;
		public readonly Vector side;

		public Direction( Vector forward, Vector top, Vector side ) {
			this.forward = forward.normalized();
			this.top     = top.normalized();
			this.side    = side.normalized();
		}

		public static Direction lookDirUp( Vector lookAt, Vector lookUp ) {
			return new Direction( lookAt, lookUp, Vector.cross( lookAt, lookUp ) );
		}

		public Direction angleAxis( Vector rotVec, double angle ) {
			return new Direction(
				Vector.angleAxis( forward, rotVec, angle ),
				Vector.angleAxis( top	 , rotVec, angle ),
				Vector.angleAxis( side   , rotVec, angle ) 
			);
		}

		public override string ToString() {
			return $"fwd : {forward}\ntop : {top}\nside: {side}";
		}
	}


	public class ManeuverNode {
		public readonly double time, radial, normal, prograde, deltaV;

		public readonly Vector burnVector;

		public readonly Orbit nextOrbit;

		public ManeuverNode( Orbit currentOrbit, double time, double radial, double normal, double prograde ) {
			this.time     = time;
			this.radial   = radial;
			this.normal   = normal;
			this.prograde = prograde;

			Direction dir_at_mn = currentOrbit.direction_at_trueAnomaly( currentOrbit.trueAnomaly_at_time(time) );
			
			burnVector = dir_at_mn.forward * prograde + dir_at_mn.top * normal + dir_at_mn.side * radial;
			deltaV = burnVector.mag();
		}

		public ManeuverNode( Orbit currentOrbit, double time, Vector scaled_burnVector ) {
			this.burnVector = scaled_burnVector;
			this.time = time;

			deltaV = burnVector.mag();

			Direction dir_at_mn = currentOrbit.direction_at_trueAnomaly( currentOrbit.trueAnomaly_at_time(time) );
			this.radial   = Vector.dot( burnVector, dir_at_mn.side );
			this.normal   =	Vector.dot( burnVector, dir_at_mn.top );
			this.prograde = Vector.dot( burnVector, dir_at_mn.forward );
		}
	}


	public class Solver {
		public static double ERROR_THRESHOLD     = 1E-5d;
		public static uint   DISCRETE_STEPS      = 5u;
		public static uint   MAX_ITERATION_STEPS = 1_000u;

		public static (double, double) newton( Func<double, double> baseFunc, Func<double, double> derivativeFunc, double startValue, double err = 1E-5f, uint maxIter = 1_000u ) {
			double currentError = Abs( 2 * err );
			double X = startValue;

			try {
				while ( currentError > err && maxIter-- >= 0 ) {
					currentError = X;
					X -= baseFunc( X ) / derivativeFunc( X );
					currentError = Abs( X - currentError );
				}
			} catch ( DivideByZeroException ) { }

			return (X, currentError);
		}

		public static (double, double) secant( Func<double, double> function, double bound0, double bound1, double err = 1E-5f, uint maxIter = 1_000u ) {
			double currentError = Abs(2 * err);

			double funcOut;
			double newBound = 0;

			try {
				while ( currentError > err && maxIter-- >= 0 ) {
					currentError = bound1;

					funcOut = function( bound1 );
					newBound = bound1 - funcOut * ( bound1 - bound0 ) / ( funcOut - function( bound0 ) );

					currentError = Abs( newBound - currentError );

					bound0 = bound1;
					bound1 = newBound;
				}
			} catch ( DivideByZeroException ) { }

			return (newBound, currentError);
		}

		public static (double, double) binary( Func<double, double> function, double left, double right, double err = 1E-5f, uint maxIter = 1_000u ) {
			double X = Abs( 2 * err );
			double midPoint = 0.0;

			while ( Abs( X ) > err && maxIter-- >= 0 ) {
				midPoint = .5f * ( left + right );
				X = function( midPoint );

				if ( X > 0 )
					right = midPoint;
				else
					left = midPoint;
			}

			return (midPoint, Abs( X ));
		}


		private static readonly Dictionary<EccAnomalyMode, Delegate> modeDictionary = new Dictionary<EccAnomalyMode, Delegate> {
			{ EccAnomalyMode.direct_bassel,	new Func<double, double, double>( bassel_eccAnomaly ) },
			{ EccAnomalyMode.direct_LIT_n5,	new Func<double, double, double>( LIT_eccAnomaly_n5 ) },
			{ EccAnomalyMode.newton_bassel,	new Func<double, double, double>( eccAnomaly_Bassel_Newton ) },
			{ EccAnomalyMode.newton_LIT_n5,	new Func<double, double, double>( eccAnomaly_LIT_Newton ) },
			{ EccAnomalyMode.newton,		new Func<double, double, double, double>( eccAnomaly_Newton ) },
			{ EccAnomalyMode.secant,		new Func<double, double, double, double, double>( eccAnomaly_Secant ) },
			{ EccAnomalyMode.binary,		new Func<double, double, double, double, double>( eccAnomaly_Binary ) },
		};
		
		public enum EccAnomalyMode {
			direct_bassel, direct_LIT_n5,
			newton_bassel, newton_LIT_n5, newton,
			secant, binary
		}

		public static Delegate get_EccAnomaly_method(EccAnomalyMode mode) {
			return modeDictionary[mode];
		}


		public static double bassel_Jn( double x, int n ) {
			double outValue = 0f;
			int kFac = 1;
			int nkFac = n + 0;

			double xHalf = .5f * x;
			double xPowN = Pow( xHalf, n );

			int signInt = -1;

			int k = 0;
			while ( k < MAX_ITERATION_STEPS ) {
				signInt *= -1;
				outValue += (signInt / (float) ( kFac * nkFac )) * xPowN * Pow( xHalf, 2 * k );
				k++;
				kFac *= k;
				nkFac *= n + k;
			}

			return outValue;
		}
		public static double bassel_eccAnomaly( double meanAnomaly, double ecc ) {
			// Bad Runtime: O( n * 5^n )
			double E = meanAnomaly; // radians

			for ( int n = 1 ; n <= DISCRETE_STEPS ; n++ ) {
				E += 2 / n * bassel_Jn( n * ecc, n ) * Sin( n * meanAnomaly );
			}

			return E;
		}

		public static double LIT_eccAnomaly_n5( double meanAnomaly, double ecc ) {
			// meanAnomaly in radians
			// Lagrange inversion theorem to solve the Eccentricity Anomaly
			// depth is only n=5, => wich results at a max Error of +-9%, for an eccentricity of 0.8

			double m = .5f * ( meanAnomaly + Constants.pi );
			double inputDiffrence = meanAnomaly - Orbit.meanAnomaly_from_eccentricityAnomaly( m, ecc );

			double C_m = Cos( m );
			double S_m = Sin( m );

			double C_2m = C_m * C_m - S_m * S_m;

			double base_denominator = ecc * C_m - 1f;
			double denominator = base_denominator * base_denominator * base_denominator;

			double c1 = 1f / ( 1f - ecc * C_m );
			double c2 = ecc * S_m / denominator;
			denominator *= base_denominator * base_denominator;

			double c3 = ecc * ecc * ( C_2m - 2 * ecc * C_m ) / denominator;
			denominator *= base_denominator * base_denominator;

			double c4 = -ecc * S_m * ( 3 * ecc * ecc * ( C_2m - 4 ) + 8 * ecc * C_m + 1 ) / denominator;

			double output = m;
			output += c1 * inputDiffrence;			// c1 / 1	= c1
			inputDiffrence *= inputDiffrence;

			output += .5f * c2 * inputDiffrence;	// c2 / 2	= .5f * c2
			inputDiffrence *= inputDiffrence;

			output += .1667f * c3 * inputDiffrence; // c3 / 6	= .1667f * c3
			inputDiffrence *= inputDiffrence;

			output += .0417f * c4 * inputDiffrence; // c4 / 24	= .0417f * c4

			return output;
		}

		public static double eccAnomaly_Bassel_Newton( double meanAnomaly, double ecc ) {
			uint newton_max_iter = MAX_ITERATION_STEPS; // the user set this to be the max steps for the newton algorithm
			MAX_ITERATION_STEPS = 5u; // sets the the bassel_Jn max steps

			(double E, double er) = newton(
				baseFunc:		( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly( E, ecc ) - meanAnomaly,
				derivativeFunc: ( double E ) => 1 - ecc*Cos( E ),
				startValue:		bassel_eccAnomaly( meanAnomaly, ecc ),
				err:			ERROR_THRESHOLD,
				maxIter:		newton_max_iter
			);
			return E;
		}

		public static double eccAnomaly_LIT_Newton( double meanAnomaly, double ecc ) {
			(double E, double er) = newton(
				baseFunc:		( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly( E, ecc ) - meanAnomaly,
				derivativeFunc: ( double E ) => 1 - ecc * Cos( E ),
				startValue:		LIT_eccAnomaly_n5( meanAnomaly, ecc ),
				err:			ERROR_THRESHOLD,
				maxIter:		MAX_ITERATION_STEPS
			);
			return E;
		}

		public static double eccAnomaly_Newton( double meanAnomaly, double ecc, double startValue=Constants.pi ) {
			(double E, double er) = newton(
				baseFunc:		( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly( E, ecc ) - meanAnomaly,
				derivativeFunc: ( double E ) => 1 - ecc * Cos( E ),
				startValue:		startValue,
				err:			ERROR_THRESHOLD,
				maxIter:		MAX_ITERATION_STEPS
			);
			return E;
		}

		public static double eccAnomaly_Secant( double meanAnomaly, double ecc, double x0=0, double x1 =Constants.tau ) {
			(double E, double er) = secant(
				function:	( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly(E, ecc) - meanAnomaly,
				bound0:		x0,
				bound1:		x1,
				err:		ERROR_THRESHOLD,
				maxIter:	MAX_ITERATION_STEPS
			);
			return E;
		}
		
		public static double eccAnomaly_Binary( double meanAnomaly, double ecc, double leftBound=0, double rightBound=Constants.tau ) {
			(double E, double er) = binary(
				function:	( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly(E, ecc) - meanAnomaly,
				left:		leftBound,
				right:		rightBound,
				err:		ERROR_THRESHOLD,
				maxIter:	MAX_ITERATION_STEPS
			);
			return E;
		}
	}
}
