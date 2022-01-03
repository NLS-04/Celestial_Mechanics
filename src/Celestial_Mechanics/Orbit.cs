using System;
using static System.Math;

using Vectors;
using Directions;

using Solver = SolverAlgorithms.Solver;

namespace Celestial_Mechanics {
	public class Constants {
		public const double G   = 6.67430E-11d;
		public const double e   = Math.E;
		public const double pi  = Math.PI;
		public const double tau = 2 * pi;

		/// <summary>Convert a degrees value into a radians value by multiplying with this constant <see cref="rad"/></summary>
		public const double rad = pi / 180;

		/// <summary>Convert a radians value into a degrees value by multiplying with this constant <see cref="deg"/></summary>
		public const double deg = 180 / pi;
	}

	/// <summary>
	/// <para>
	/// Frame Of Refrences:<br />
	///		=>	ECI		('Earth'-Centered-Inertial): is a inertial Refrence Frame, on which all orbit parameters are defined<br />
	///		=>	ECEF	('Earth'-Centered-'Eaerth'-Fixed): is a non-inertial Refrence Frame, which is usefull to describe positions relative to a Body's surface<br />
	/// </para>
	/// </summary>
	public class Body {
		// mu = G*M
		/// <summary>Standard gravitational parameter 'mu' of an <see cref="Body"/> in [m^3/s^2]</summary>
		public readonly double mu;
		/// <summary>Mass of an <see cref="Body"/> in [kg]</summary>
		public readonly double mass;
		/// <summary>Average radial distance of the surface of an <see cref="Body"/> in [m]</summary>
		public readonly double radius;


		public  Orbit? orbit { get => __orbit; }
		private Orbit? __orbit;

		public  double? sphereOfInfluence { get => __sphereOfInfluence; }
		private double? __sphereOfInfluence;


		public readonly double rotationPeriod;	// time in [seconds] that it takes to rotate this body by 360°
		public readonly double rotationSpeed;	// rate of change of angle in [radians] per time in [seconds] => [radians/second]
		public readonly Vector rotationVector;	// normalized Vector pointing in the direction of rotation

		/// <summary>This <see cref="Direction"/> is defined to be this Body's ECI Frame of refrence at <see cref="Orbit.TimeMode.UT"/> <see langword="0.0d"/>, thus all derived orbit parameters are based on this (inertial) ECI Frame </summary>
		public readonly Direction ECI_Orientation = Direction.defaultDirection;


		/* Example of how to Convert Vectors between diffrent Bodies (Coordinate Systems):
		 *		!! LOOKUP THE PNG (./src/coordinateConvertion.png) FOR VISUALISATION AND BETTER UNDERSTANDING !!
		 * 
		 *	-----------------------------------------------------------------------------------------------------------------------------------------------------------
		 *	
		 *	1. Convert a Satelite position from Earths ECI to Moons ECI refrence frame. (assuming Moon is a (orbiting) sub Body of Earth)
		 *		Concept:
		 *			( [Earth ECI, Moon as Origin] pos_Satelite ) = ( [Earth ECI, Earth as Origin] pos_Satelite ) - ( [Earth ECI, Earth as Origin] pos_Moon )
		 *			( [Moon ECI, Moon as Origin] pos_Satelite ) = Moon.WorldToLocal( ([Earth ECI, Moon as Origin] pos_Satelite) )
		 *		
		 *		Pseudo-ish Code:
		 *			Vector pos_Satelite = Earth.getPositionOf( Satelite, TimeNow ) - Earth.getPositionOf( Moon, TimeNow );
		 *			return Moon.WorldToLocal( pos_Satelite );
		 *		
		 *	-----------------------------------------------------------------------------------------------------------------------------------------------------------
		 *	
		 *	2. Convert a Satelite position from Moons ECI to Earths ECI refrence frame. (assuming Moon is a (orbiting) sub Body of Earth)
		 *		Concept:
		 *			( [Earth ECI] pos_Satelite ) = Earth.worldToLocal( ([Moon ECI] pos_Satelite) )
		 *			( [Earth ECI] pos_Satelite ) = ( [Earth ECI] pos_Satelite ) + ( [Earth ECI] pos_Moon )
		 *		
		 *		Pseudo-ish Code:
		 *			Vector pos_Satelite = Earth.worldToLocal( pos_Satelite );
		 *			return pos_Satelite + Earth.getPositionOf( Moon, TimeNow_UT );
		 *			
		 *	-----------------------------------------------------------------------------------------------------------------------------------------------------------
		 */


		/// <summary>Constructor for an (large) orbitable Body</summary>
		/// <param name="mu">Standard gravitational parameter 'mu' of <see langword="this"/> <see cref="Body"/> in [m^3/s^2]</param>
		/// <param name="radius">Average radial distance of the surface of <see langword="this"/> <see cref="Body"/> in [m]</param>
		/// <param name="orient">A struct to orientate a coordinate system of <see langword="this"/> <see cref="Body"/>; It is defined as being the <see cref="Direction"/> of <see langword="this"/> <see cref="Body"/> at Uinversal Time (<see cref="Orbit.TimeMode.UT"/>) <see langword="0.0d"/>; If its <see langword="null"/> the <see cref="Body.orientation"/> is set to the <see cref="defaultDirection"/> </param>
		/// <param name="angularVelocity">A <see cref="Vector"/> pointing in the Direction of the Axis of rotation; The magnitude of the <see cref="Vector"/> is the speed of the rotation in [radians / second]; If its <see langword="null"/> <see cref="Vector.zero"/> is being assigned</param>
		public Body( float mu, double radius, Vector? angularVelocity=null, Direction? orient=null ) {
			this.mu     = mu;
			this.mass   = mu / Constants.G;
			this.radius = radius;
			
			this.rotationVector = angularVelocity ?? Vector.zero;
			this.ECI_Orientation = orient ?? Direction.defaultDirection;

			this.rotationSpeed = this.rotationVector.mag();
			this.rotationPeriod = Constants.tau / this.rotationSpeed;
		}

		/// <summary>Constructor for an (large) orbitable Body</summary>
		/// <param name="radius">Average radial distance of the surface of <see langword="this"/> <see cref="Body"/> in [m]</param>
		/// <param name="mass">Mass of <see langword="this"/> <see cref="Body"/> in [kg]</param>
		/// <param name="orient">A struct to orientate a coordinate system of <see langword="this"/> <see cref="Body"/>; It is defined as being the <see cref="Direction"/> of <see langword="this"/> <see cref="Body"/> at Uinversal Time (<see cref="Orbit.TimeMode.UT"/>) <see langword="0.0d"/>; If its <see langword="null"/> the <see cref="Body.orientation"/> is set to the <see cref="defaultDirection"/> </param>
		/// <param name="angularVelocity">A <see cref="Vector"/> pointing in the Direction of the Axis of rotation; The magnitude of the <see cref="Vector"/> is the speed of the rotation in [radians / second]; If its <see langword="null"/> <see cref="Vector.zero"/> is being assigned</param>
		public Body( double radius, float mass, Vector? angularVelocity=null, Direction? orient=null ) {
			this.mass   = mass;
			this.mu     = mass * Constants.G;
			this.radius = radius;
			
			this.rotationVector = angularVelocity ?? Vector.up;
			this.ECI_Orientation = orient ?? Direction.defaultDirection;

			this.rotationSpeed = this.rotationVector.mag();
			this.rotationPeriod = Constants.tau / this.rotationSpeed;
		}


		public readonly static Body KERBIN = new Body( 3.5316E12f,      600_000d  , Vector.up * Constants.tau / 21_549.425d );
		public readonly static Body EARTH  = new Body( 3.986004418E14f, 6_371_000d, Vector.up * Constants.tau / 86_164.0d   );
		public readonly static Body MARS   = new Body( 4.282837E13f,	3_389_500d, Vector.up * Constants.tau / 88_642.0d   );


		public override bool Equals( object obj ) {
			if ( !(obj is Body) )
				return false;
			Body bod = (Body) obj;
			return mu == bod.mu && radius == bod.radius && ( ECI_Orientation.Equals(bod.ECI_Orientation) );
		}
		public override int GetHashCode() => -100;

		public static bool operator ==( Body a, Body b ) => a.Equals( b );
		public static bool operator !=( Body a, Body b ) =>!a.Equals( b );

		public void setOrbit( Orbit orb ) {
			__orbit = orb;
			__sphereOfInfluence = orbit.semiMajorAxis * Pow( mass / orbit.body.mass, 2d / 5d );
		}


		public Direction getDirectionAtTime( double UT ) {
			// the modulo is not actually necessary, but clarifies that the angle should be between 0 and 2*pi
			return ECI_Orientation.angleAxis( rotationVector, ( rotationSpeed * UT ) % Constants.tau );
		}

		public Vector getParentPosition( double UT ) {
			if ( orbit is null )
				return Vector.zero;
			return orbit.position_at_time( UT, Orbit.TimeMode.UT );
		}

		public Vector getWorldPosition( double UT ) {
			return convert_LocalToWorld( UT, Vector.zero );
		}

		public Vector convert_LocalToWorld( double UT, Vector p ) {
			if ( orbit is null ) // we are the center of the univers :-) so the vector is allready in the correct Frame
				return p;

			// convert this Vector definded in this Body's coordinate system to the system of this Body's parent Body
			Vector convertedPos = orbit.body.ECI_Orientation.worldToLocal( p );
			
			// add the position of this Body in the frame of this Body's parent Body
			convertedPos += orbit.position_at_time( UT, Orbit.TimeMode.UT );

			return orbit.body.convert_LocalToWorld( UT, convertedPos );
		}

		public Vector convert_WorldToLocal( double UT, Vector p ) {
			Vector offsetedVector = p - getWorldPosition( UT );

			return ECI_Orientation.worldToLocal( offsetedVector );
		}

		/// <summary>Find the <see cref="Body"/> that connects both <see cref="Body"/>s <paramref name="a"/> and <paramref name="b"/> via the shortest Path</summary>
		public static Body get_same_parent_body( Body a, Body b ) {
			(Orbit tmp_a, Orbit tmp_b) = (a.orbit, b.orbit);

			while ( tmp_a != null ) {
				while ( tmp_b != null ) {
					if ( tmp_a.body == tmp_b.body )
						return tmp_a.body;

					tmp_b = tmp_b.body.orbit;
				}

				tmp_a = tmp_a.body.orbit;
			}

			return a;
		}
	}

	/// <summary>
	///	<para>
	///	-> LAN (longitude of Ascending Node) measured to X-Axis of the Body's initial ECI-Frame<br />
	/// -> Inclination measured to XY-Plane of the Body's initial ECI-Frame<br />
	/// -> Body's North == Z-Axis of the Body's initial ECI-Frame
	/// </para>
	/// 
	/// <remarks>
	/// All (readable) angle Values are in [DEGREES] but are internally converted and calculated as [RADIANS]
	/// </remarks>
	/// 
	/// </summary>
	public class Orbit {

		public Body body;

		/// <summary>The Vector is normalized, and its actual magnitude is separately stored at an associativ variable</summary>
		public readonly Vector eccentricity_Vector, angularMomentum_Vector, nodeLine_Vector;

		public readonly double apoapsis, periapsis, semiMajorAxis, semiMinorAxis;
		public readonly double inclination, eccentricity, angularMomentum, longitudeOfAscendingNode, argumentOfPeriapsis;
		public readonly double period, epoch, meanAnomaly_At_Epoch, semi_latus_rectum;

		public ManeuverNode? maneuverNode { get; set; }

		/// <summary>
		/// <see cref="Direction.forward"/> = prograde Direction<br/>
		/// <see cref="Direction.side"/> = radial (out) Direction<br/>
		/// <see cref="Direction.top"/> = normal Direction
		/// </summary>
		private readonly Direction dir_at_periapsis;

		#region constructors
		/// <summary><see cref="Orbit"/> constructor for general Orbits, orientated in 3D</summary>
		/// <param name="inc">Inclination in [degrees]</param>
		/// <param name="ecc">Eccentricity of Orbit, ecc >= 0</param>
		/// <param name="sma">Semi-Major-Axis of Orbit in [m]</param>
		/// <param name="lan">Longitude of Ascending Node in [degrees]</param>
		/// <param name="argOfPer">Argument of Periapsis in [degrees]</param>
		/// <param name="mEp">Mean Anomaly in [Degrees] at Epoch</param>
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


			nodeLine_Vector        = Vector.angleAxis( body.ECI_Orientation.forward, body.ECI_Orientation.top, longitudeOfAscendingNode * Constants.rad ).normalized();
			angularMomentum_Vector = Vector.angleAxis( body.ECI_Orientation.top, nodeLine_Vector, inclination * Constants.rad ).normalized();
			eccentricity_Vector    = Vector.angleAxis( nodeLine_Vector, angularMomentum_Vector, argOfPer * Constants.rad ).normalized();

			meanAnomaly_At_Epoch = mEp;

			dir_at_periapsis = new Direction( -Vector.cross( eccentricity_Vector, angularMomentum_Vector ), eccentricity_Vector, angularMomentum_Vector );
		}

		/// <summary><see cref="Orbit"/> constructor for general Orbits, orientated in 3D, based on an initial position and velocity</summary>
		/// <param name="pos">initial position of satellite in [m] in Body's [ECI]-Frame</param>
		/// <param name="vel">initial velocity of satellite in [m/s] in Body's [ECI]-Frame</param>
		/// <param name="tEpoch">Arbitrary timestamp to calculate time dependend parameters of this orbit (e.g. position, velocity, flight angle etc.)</param>
		/// <param name="body"><see cref="Body"/> to be orbited around</param>
		public Orbit( Vector pos, Vector vel, double tEpoch, Body body ) {
			epoch     = tEpoch;
			this.body = body;

			angularMomentum_Vector = Vector.cross( pos, vel ).normalized();
			angularMomentum        = angularMomentum_Vector.mag();

			inclination = Vector.vang( angularMomentum_Vector, Vector.up ) * Constants.deg;

			nodeLine_Vector = Vector.cross( Vector.up, angularMomentum_Vector ).normalized();

			longitudeOfAscendingNode = Atan2( nodeLine_Vector.y, nodeLine_Vector.x ) * Constants.deg;

			eccentricity_Vector = ( angularMomentum / body.mu ) * ( Vector.cross( vel, angularMomentum_Vector ) - body.mu * pos.normalized() );
			eccentricity        = eccentricity_Vector.mag();

			argumentOfPeriapsis = Vector.vang( nodeLine_Vector, eccentricity_Vector ) * Constants.deg;
			argumentOfPeriapsis = eccentricity_Vector.z < 0 ? 360 - argumentOfPeriapsis : argumentOfPeriapsis;

			apoapsis  = altitude( 180 );
			periapsis = altitude(   0 );

			semiMajorAxis = .5d * ( apoapsis + periapsis );
			semiMinorAxis = calculate_semiMinorAxis( semiMajorAxis, eccentricity );
			period        = calculate_Period( semiMajorAxis, body.mu );

			semi_latus_rectum    = angularMomentum * angularMomentum / body.mu;
			eccentricity_Vector  = eccentricity_Vector.normalized();
			meanAnomaly_At_Epoch = meanAnomaly( Vector.vang( body.ECI_Orientation.forward, eccentricity_Vector ) * Constants.deg );

			dir_at_periapsis = new Direction( -Vector.cross( eccentricity_Vector, angularMomentum_Vector ), eccentricity_Vector, angularMomentum_Vector );
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
		public override int GetHashCode() => -500;
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
			return period * meanAnomaly * Constants.rad / Constants.tau;
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

		#region Two Orbit Methods
		public struct PointsOfInterest {
			/// <summary>I, J being the angles in [degrees] of interest</summary>
			public double? I, J;
			public bool hasValues { get => I.HasValue && J.HasValue; }
		}

		public static PointsOfInterest solveTheta( double a, double b, double c ) {
			// solve the angle theta in [degrees] for the equation of the form:
			// a * cos(theta) + b * sin(theta) = c

			if ( a == 0d )
				return new PointsOfInterest { };

			double phi = Atan( b / a );
			double cosContent = c / a * Cos(phi);

			if ( cosContent < -1d || cosContent > 1d )
				return new PointsOfInterest { };

			cosContent = Acos( cosContent );
			
			return new PointsOfInterest { 
				I = ( phi + cosContent ) * Constants.deg, 
				J = ( phi - cosContent ) * Constants.deg 
			};
		}

		/// <summary>
		/// Calculate the Points there the <see cref="Orbit"/> <paramref name="A"/> has an intersection with <see cref="Orbit"/> <paramref name="B"/><br/>
		/// implying that:<br/>
		///		- A and B are Orbiting the same <see cref="Body"/><br/>
		///		- A and B are in the same (orbital) plane<br/>
		/// </summary>
		/// <returns>
		/// The True Anomalies as <see cref="PointsOfInterest"/> in respect to the <see cref="Orbit"/> <paramref name="A"/><br/>
		///	There is NOT allways a Solution (<see cref="PointsOfInterest.hasValues"/> == <see langword="true"/>)
		/// </returns>
		public static PointsOfInterest getIntersectionPoints( Orbit A, Orbit B ) {
			if ( A.body != B.body )
				throw new ArgumentException( "Orbit A and B MUST have the same Body which they orbit" );

			double eta = Constants.rad * (B.argumentOfPeriapsis - A.argumentOfPeriapsis);

			PointsOfInterest poi = solveTheta(
				 A.eccentricity * B.angularMomentum * B.angularMomentum - B.eccentricity * A.angularMomentum * A.angularMomentum * Cos( eta ),
				-B.eccentricity * A.angularMomentum * A.angularMomentum * Sin( eta ),
				 A.angularMomentum * A.angularMomentum - B.angularMomentum * B.angularMomentum
			);

			poi.I = ( poi.I + 360 ) % 360;
			poi.J = ( poi.J + 360 ) % 360;

			return poi;
		}

		/// <summary>
		/// Calculate the Points there the <see cref="Orbit"/> <paramref name="A"/> is tangential to <see cref="Orbit"/> <paramref name="B"/><br/>
		/// implying that:<br/>
		///		- A and B are Orbiting the same <see cref="Body"/>
		/// </summary>
		/// <returns>
		/// The True Anomalies as <see cref="PointsOfInterest"/> in respect to the <see cref="Orbit"/> <paramref name="A"/><br/>
		///	There IS allways a Solution (<see cref="PointsOfInterest.hasValues"/> == <see langword="true"/>)
		/// </returns>
		public static PointsOfInterest getTangentialPoints( Orbit A, Orbit B ) {
			if ( A.body != B.body )
				throw new ArgumentException( "Orbit A and B MUST have the same Body which they orbit" );

			double eta = Constants.rad * (B.argumentOfPeriapsis - A.argumentOfPeriapsis);

			PointsOfInterest poi = solveTheta(
				B.eccentricity * Sin(eta),
				A.eccentricity - B.eccentricity * Cos(eta),
				A.eccentricity * B.eccentricity * Sin( eta )
			);

			poi.I = ( poi.I + 360 ) % 360;
			poi.J = ( poi.J + 360 ) % 360;

			return poi;
		}


		public static double getDistance_at_Time( Orbit A, Orbit B, double UT ) {
			(Orbit orbitA, Orbit orbitB) = (A.get_Orbit_at_time( UT ), B.get_Orbit_at_time( UT ));

			Vector posA, posB;
			if ( orbitA.body == orbitB.body ) {
				posA = orbitA.position_at_time( UT, TimeMode.UT );
				posB = orbitB.position_at_time( UT, TimeMode.UT );
			} else {
				posA = A.position_world_at_time( UT );
				posB = B.position_world_at_time( UT );
			}

			return ( posA - posB ).mag();
		}

		public static (double time, Vector pos) get_closest_Approach( Orbit A, Orbit B, (double start, double range) timeBoundary/*, double maxDistance=10_000d*/ ) {
			Solver.Result res = Solver.section( (double UT) => Orbit.getDistance_at_Time( A, B, UT ), timeBoundary.start, timeBoundary.range, 3 );
			return ( (double) res.input, A.position_world_at_time( (double) res.input ) );
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

		/// <summary>calculating the <see cref="GeoPosition"/> of an Object on this initial <see cref="Orbit"/> at a given Time</summary>
		/// <param name="UT">Time in [seconds] in [<see cref="TimeMode.UT"/>]</param>
		/// <returns>The <see cref="GeoPosition"/> at this UT</returns>
		public GeoPosition getGeoPosition( double UT ) {
			return new GeoPosition( position_at_time( UT, TimeMode.UT ), body );
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

		/// <summary>Calculates the time it takes to transvers between two True Anomalies for an Satelite on this Orbit</summary>
		/// <param name="tA_from">Angle in [degrees] of an position to the periapsis</param>
		/// <param name="tA_to">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>the delta Time in [seconds]</returns>
		public double get_deltaTime_of_trueAnomalies( double tA_from, double tA_to ) {
			double dT = trueAnomaly_to_time( tA_to ) - trueAnomaly_to_time( tA_from );
			return ( dT + period ) % period;
		}

		
		/// <summary>Calculate the Mean Anomaly at a given time</summary>
		/// <param name="time">Time in [seconds]</param>
		/// <param name="TM">The Time mode controls how the time value should be interpreted</param>
		/// <returns>The Mean Anomaly in [degrees]</returns>
		public double meanAnomaly_at_time( double time, TimeMode TM=TimeMode.UT ) {
			double Me = 0d;
			Me = TM switch 
			{
				TimeMode.UT => Constants.tau * ( time - epoch ) / period + meanAnomaly_At_Epoch,
				TimeMode.PeriapsRelativ => Constants.tau * time / period,
				_ => throw new ArgumentException("Use Correct TimeMode, e.g. UT, PeriapsRelativ", "TM")
			};

			return ( Me % Constants.tau ) * Constants.deg;
		}

		/// <summary>
		/// Calculate the True Anomaly at a given time<br/>
		/// The solution is non trivial, thus this value gets calculated via an iterative Algorithm</summary>
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

		/// <summary>Calculating the <see cref="Orbit"/> an Object on this initial <see cref="Orbit"/> will end up, after performing all up coming (up to that UT) Maneuver Events </summary>
		/// <param name="UT">Time in [seconds] in [<see cref="TimeMode.UT"/>]</param>
		/// <returns>The <see cref="Orbit"/> at this UT</returns>
		public Orbit get_Orbit_at_time( double UT ) {
			Orbit finalOrbit = this;

			while ( true ) {
				if ( finalOrbit.maneuverNode == null )
					break;

				if ( UT < finalOrbit.maneuverNode.time )
					break;

				finalOrbit = finalOrbit.maneuverNode.nextOrbit;
			}

			return finalOrbit;
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
		/// <returns>The Position 3d-normalized-Vector in (this body's) [ECI]-Frame</returns>
		public Vector normalized_position( double trueAnomaly ) {
			// r_hat or u_hat, position direction
			return Vector.angleAxis( eccentricity_Vector, angularMomentum_Vector, trueAnomaly ).normalized();
		}


		/// <summary>Calculate the Position Vector at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The Position 3d-Vector in [meters] in (this body's) [ECI]-Frame</returns>
		public Vector position_at_trueAnomaly( double trueAnomaly ) {
			return altitude( trueAnomaly ) * normalized_position( trueAnomaly );
		}

		/// <summary>Calculate the Velocity Vector at a given True Anomaly</summary>
		/// <param name="trueAnomaly">Angle in [degrees] of an position to the periapsis</param>
		/// <returns>The Velocity 3d-Vector [meters/second] in (this body's) [ECI]-Frame</returns>
		public Vector velocity_at_trueAnomaly( double trueAnomaly ) {
			Direction dir = direction_at_trueAnomaly( trueAnomaly );
			Vector radial_out_vel = radial_velocity( trueAnomaly ) * dir.side;
			Vector tangential_vel = tangential_velocity( trueAnomaly ) * dir.forward;
			return radial_out_vel + tangential_vel;
		}


		/// <summary>
		/// Calculate the Position Vector at a given Time<br/>
		/// Is the <paramref name="TM"/> =<see cref="TimeMode.UT"/> then all up coming (up to that UT) Maneuver Events are being regared
		/// </summary>
		/// <param name="time">Time in [seconds]</param>
		/// <param name="TM">The Time mode controls how the time value should be interpreted</param>
		/// <returns>The Position 3d-Vector in [meters] in (this body's) [ECI]-Frame</returns>
		public Vector position_at_time( double time, TimeMode TM = TimeMode.UT ) {
			switch ( TM ) {
				case TimeMode.UT:
					Orbit finalOrbit = get_Orbit_at_time( time );
					return finalOrbit.position_at_trueAnomaly( finalOrbit.trueAnomaly_at_time( time, TimeMode.UT ) );
			
				case TimeMode.PeriapsRelativ:
					return position_at_trueAnomaly( trueAnomaly_at_time( time, TimeMode.PeriapsRelativ ) );
				
				default:
					return null;
			}
		}

		/// <summary>
		/// Calculate the Velocity Vector at a given Time<br/>
		/// Is the <paramref name="TM"/> =<see cref="TimeMode.UT"/> then all up coming (up to that UT) Maneuver Events are being regared
		/// </summary>
		/// <param name="time">Time in [seconds]</param>
		/// <param name="TM">The Time mode controls how the time value should be interpreted</param>
		/// <returns>The Velocity 3d-Vector [meters/second] in (this body's) [ECI]-Frame</returns>
		public Vector velocity_at_time( double time, TimeMode TM=TimeMode.UT ) {
			switch ( TM ) {
				case TimeMode.UT:
					Orbit finalOrbit = get_Orbit_at_time( time );
					return finalOrbit.velocity_at_trueAnomaly( finalOrbit.trueAnomaly_at_time( time, TimeMode.UT ) );

				case TimeMode.PeriapsRelativ:
					return velocity_at_trueAnomaly( trueAnomaly_at_time( time, TimeMode.PeriapsRelativ ) );

				default:
					return null;
			}
		}


		/// <summary>
		/// Calculate the World-Position-Vector at a given Time<br/>
		/// Regarding all up coming (up to that <paramref name="UT"/>) Maneuver Events
		/// </summary>
		/// <param name="UT">Time in [seconds] in [<see cref="TimeMode.UT"/>]</param>
		/// <returns>The World-Position 3d-Vector in [meters] in [World]-Frame</returns>
		public Vector position_world_at_time( double UT ) {
			Orbit orb = get_Orbit_at_time( UT );
			return orb.body.convert_LocalToWorld( UT, orb.position_at_time( UT, TimeMode.UT ) );
		}

		/// <summary>
		/// Calculate the World-Velocity-Vector at a given Time<br/>
		/// Regarding all up coming (up to that <paramref name="UT"/>) Maneuver Events
		/// </summary>
		/// <param name="UT">Time in [seconds] in [<see cref="TimeMode.UT"/>]</param>
		/// <returns>The World-Velocity 3d-Vector in [meters] in [World]-Frame</returns>
		public Vector velocity_world_at_time( double UT ) {
			Orbit orb = get_Orbit_at_time( UT );
			return orb.body.convert_LocalToWorld( UT, orb.velocity_at_time( UT, TimeMode.UT ) );
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

	public struct GeoPosition {
		// (latitude  "=") declination     = 'Vertical location measure'
		// (longitude "=") right ascension = 'Horizontal location measure'

		public readonly Body body;
		/// <summary>Time independent Geoposition angle measurement in [degrees], in body's [ECI]-Frame</summary>
		public readonly double declination, rightAscension;

		public readonly Vector ECI_vector;

		public GeoPosition( double dec, double rightAsc, Body bod ) {
			body = bod;
			declination = dec * Constants.deg;
			rightAscension = rightAsc * Constants.deg;

			Vector decVec = Vector.angleAxis( body.ECI_Orientation.forward, body.rotationVector, dec );
			ECI_vector = Vector.angleAxis( decVec, body.ECI_Orientation.top, rightAsc );
		}

		/// <param name="pos">The position Vector MUST be in bods's [ECI]-Frame</param>
		public GeoPosition( Vector pos, Body bod ) {
			body = bod;

			ECI_vector = pos.normalized();
			(declination, rightAscension) = GeoPosition.declination_RightAscension( ECI_vector );
		}

		/// <summary>Latitude, Longitude angles at a given Time, e.g. Converts the declination and right Ascension ECI angles into ECEF angles</summary>
		/// <param name="UT"><see cref="Orbit.TimeMode.UT"/></param>
		/// <returns>(Latitude, Longitude) in [degrees], in body's [ECEF]-Frame</returns>
		public (double, double) getLatLngAtTime( double UT ) {
			Vector ECEF = body.getDirectionAtTime( UT ).worldToLocal( ECI_vector );

			return GeoPosition.declination_RightAscension( ECEF );
		}

		private static (double, double) declination_RightAscension( Vector pos ) {
			double declination = Asin( pos.y ) * Constants.deg;
			double rightAscension = Acos( pos.x / Cos( declination ) ) * Constants.deg;
			return (declination, rightAscension);
		}
	}

	public class ManeuverNode {
		public readonly double time, radial, normal, prograde, deltaV;

		public readonly Vector burnVector;

		public readonly Orbit nextOrbit;

		public ManeuverNode( Orbit currentOrbit, double UT, double radial, double normal, double prograde ) {
			this.time     = UT;
			this.radial   = radial;
			this.normal   = normal;
			this.prograde = prograde;

			Direction dir_at_mn = currentOrbit.direction_at_trueAnomaly( currentOrbit.trueAnomaly_at_time(UT) );
			
			burnVector = dir_at_mn.forward * prograde + dir_at_mn.top * normal + dir_at_mn.side * radial;
			deltaV = burnVector.mag();

			nextOrbit = new Orbit( currentOrbit.position_at_time( UT, Orbit.TimeMode.UT ), burnVector + currentOrbit.velocity_at_time( UT, Orbit.TimeMode.UT ), UT, currentOrbit.body );
		}

		public ManeuverNode( Orbit currentOrbit, double UT, Vector scaled_burnVector ) {
			this.burnVector = scaled_burnVector;
			this.time = UT;

			deltaV = burnVector.mag();

			Direction dir_at_mn = currentOrbit.direction_at_trueAnomaly( currentOrbit.trueAnomaly_at_time(UT) );
			(prograde, radial, normal) = dir_at_mn.worldToLocal( burnVector );

			nextOrbit = new Orbit( currentOrbit.position_at_time( UT, Orbit.TimeMode.UT ), burnVector + currentOrbit.velocity_at_time( UT, Orbit.TimeMode.UT ), UT, currentOrbit.body );
		}
	}
}