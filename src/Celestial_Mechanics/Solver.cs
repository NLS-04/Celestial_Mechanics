using static System.Math;
using System;
using System.Collections.Generic;
using Celestial_Mechanics;

namespace SolverAlgorithms {
	public static class Solver {
		public static double ERROR_THRESHOLD     = 1E-5d;
		public static uint   DISCRETE_STEPS      = 5u;
		public static uint   MAX_ITERATION_STEPS = 1_000u;

		public struct Result {
			public readonly double  absoluteError;
			public readonly double  value;
			public readonly double? input;

			public Result( double err, double val ) {
				(absoluteError, value, input) = (err, val, null);
			}

			public Result( double err, double val, double inpu ) {
				(absoluteError, value, input) = (err, val, inpu);
			}

			public override string ToString() => $"val:\t{value}\nerror:\t{absoluteError}" + ( input is null ? "" : $"\ninput:\t{input}" );
		}

		public static Result newton( Func<double, double> baseFunc, Func<double, double> derivativeFunc, double startValue, double err = 1E-5d, uint maxIter = 1_000u ) {
			double currentError = Abs( 2 * err );
			double X = startValue;

			try {
				while ( currentError > err && maxIter-- >= 0 ) {
					currentError = X;
					X -= baseFunc( X ) / derivativeFunc( X );
					currentError = Abs( X - currentError );
				}
			} catch ( DivideByZeroException ) { }

			return new Result( X, currentError );
		}

		public static Result secant( Func<double, double> function, double bound0, double bound1, double err = 1E-5d, uint maxIter = 1_000u ) {
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

			return new Result( newBound, currentError );
		}

		public static Result binary( Func<double, double> function, double left, double right, double err = 1E-5d, uint maxIter = 1_000u ) {
			double X = Abs( 2 * err );
			double midPoint = 0.0;

			while ( Abs( X ) > err && maxIter-- >= 0 ) {
				midPoint = .5d * ( left + right );
				X = function( midPoint );

				if ( X > 0 )
					right = midPoint;
				else
					left = midPoint;
			}

			return new Result( midPoint, Abs( X ) );
		}

		public static Result section( Func<double, double> function, double start, double range, int steps, bool searchingForMin = true, double err = 1E-5d ) {
			double bestInput, bestFuncValue, input, funcVal, stepSize;
			bestInput = start;
			bestFuncValue = searchingForMin ? double.MaxValue : double.MinValue;

			/*	since the range gets readjusted after every iteration to:
			 *		new_range = old_range / steps
			 *	=>	range_of_nth_steps = inital_range / pow( steps, n )
			 *	
			 *	there exists a n (element of the Naturals) so that:
			 *		err > inital_range / pow( steps, n )
			 *	with err (element of the Real) being the Error and being greater than the range of the nth iteration
			 *	
			 *	thus,
			 *	=>	pow( steps, n ) > inital_range / err
			 *	=>	n * ln( steps ) > ln( inital_range / err )
			 *	=>	n > ln( inital_range / err ) / ln( steps )
			 *	
			 *	=>	n >= ceil( ln( inital_range / err ) / ln( steps ) )
			 *	
			 *	Runtime is therefor roughtly: O( n / ln(n) )
			 */
			int n = (int) Ceiling( Log( range / err ) / Log( steps ) );

			for ( int i = 0 ; i < n ; i++ ) {
				stepSize = range / steps;

				for ( int step = 0 ; step <= steps ; step++ ) {
					input = start + step * stepSize;
					funcVal = function( input );

					if ( ( searchingForMin && funcVal < bestFuncValue ) || ( !searchingForMin && funcVal > bestFuncValue ) ) {
						(bestInput, bestFuncValue) = (input, funcVal);
					}
				}

				// the adjusted search region is now centred around the best input and has now a range of the previous stepSize
				start = bestInput - .5d * stepSize;
				range = stepSize;
			}

			// ...*steps, since range is currently equal to the previous stepsize, and thus must be rescaled to be the prvious range
			return new Result( Abs( range * steps ), bestFuncValue, bestInput );
		}


		private static readonly Dictionary<EccAnomalyMode, Delegate> modeDictionary = new Dictionary<EccAnomalyMode, Delegate> {
			{ EccAnomalyMode.direct_bassel, new Func<double, double, double>( bassel_eccAnomaly ) },
			{ EccAnomalyMode.direct_LIT_n5, new Func<double, double, double>( LIT_eccAnomaly_n5 ) },
			{ EccAnomalyMode.newton_bassel, new Func<double, double, double>( eccAnomaly_Bassel_Newton ) },
			{ EccAnomalyMode.newton_LIT_n5, new Func<double, double, double>( eccAnomaly_LIT_Newton ) },
			{ EccAnomalyMode.newton,        new Func<double, double, double, double>( eccAnomaly_Newton ) },
			{ EccAnomalyMode.secant,        new Func<double, double, double, double, double>( eccAnomaly_Secant ) },
			{ EccAnomalyMode.binary,        new Func<double, double, double, double, double>( eccAnomaly_Binary ) },
		};

		public enum EccAnomalyMode {
			direct_bassel, direct_LIT_n5,
			newton_bassel, newton_LIT_n5, 
			newton, secant, binary
		}

		public static Delegate get_EccAnomaly_method( EccAnomalyMode mode ) {
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
				outValue += ( signInt / (float) ( kFac * nkFac ) ) * xPowN * Pow( xHalf, 2 * k );
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

			double m = .5d * ( meanAnomaly + Constants.pi );
			double inputDiffrence = meanAnomaly - Orbit.meanAnomaly_from_eccentricityAnomaly( m, ecc );

			double C_m = Cos( m );
			double S_m = Sin( m );

			double C_2m = C_m * C_m - S_m * S_m;

			double base_denominator = ecc * C_m - 1d;
			double denominator = base_denominator * base_denominator * base_denominator;

			double c1 = 1d / ( 1d - ecc * C_m );
			double c2 = ecc * S_m / denominator;
			denominator *= base_denominator * base_denominator;

			double c3 = ecc * ecc * ( C_2m - 2 * ecc * C_m ) / denominator;
			denominator *= base_denominator * base_denominator;

			double c4 = -ecc * S_m * ( 3 * ecc * ecc * ( C_2m - 4 ) + 8 * ecc * C_m + 1 ) / denominator;

			double output = m;
			output += c1 * inputDiffrence;          // c1 / 1	= c1
			inputDiffrence *= inputDiffrence;

			output += .5d * c2 * inputDiffrence;    // c2 / 2	= .5d * c2
			inputDiffrence *= inputDiffrence;

			output += .1667d * c3 * inputDiffrence; // c3 / 6	= .1667d * c3
			inputDiffrence *= inputDiffrence;

			output += .0417d * c4 * inputDiffrence; // c4 / 24	= .0417d * c4

			return output;
		}

		public static double eccAnomaly_Bassel_Newton( double meanAnomaly, double ecc ) {
			uint newton_max_iter = MAX_ITERATION_STEPS; // the user set this to be the max steps for the newton algorithm
			MAX_ITERATION_STEPS = 5u; // sets the the bassel_Jn max steps

			Result re = newton(
				baseFunc:       ( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly( E, ecc ) - meanAnomaly,
				derivativeFunc: ( double E ) => 1 - ecc*Cos( E ),
				startValue:     bassel_eccAnomaly( meanAnomaly, ecc ),
				err:            ERROR_THRESHOLD,
				maxIter:        newton_max_iter
			);
			return re.value;
		}

		public static double eccAnomaly_LIT_Newton( double meanAnomaly, double ecc ) {
			Result re = newton(
				baseFunc:       ( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly( E, ecc ) - meanAnomaly,
				derivativeFunc: ( double E ) => 1 - ecc * Cos( E ),
				startValue:     LIT_eccAnomaly_n5( meanAnomaly, ecc ),
				err:            ERROR_THRESHOLD,
				maxIter:        MAX_ITERATION_STEPS
			);
			return re.value;
		}

		public static double eccAnomaly_Newton( double meanAnomaly, double ecc, double startValue = Constants.pi ) {
			Result re = newton(
				baseFunc:       ( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly( E, ecc ) - meanAnomaly,
				derivativeFunc: ( double E ) => 1 - ecc * Cos( E ),
				startValue:     startValue,
				err:            ERROR_THRESHOLD,
				maxIter:        MAX_ITERATION_STEPS
			);
			return re.value;
		}

		public static double eccAnomaly_Secant( double meanAnomaly, double ecc, double x0 = 0, double x1 = Constants.tau ) {
			Result re = secant(
				function:   ( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly(E, ecc) - meanAnomaly,
				bound0:     x0,
				bound1:     x1,
				err:        ERROR_THRESHOLD,
				maxIter:    MAX_ITERATION_STEPS
			);
			return re.value;
		}

		public static double eccAnomaly_Binary( double meanAnomaly, double ecc, double leftBound = 0, double rightBound = Constants.tau ) {
			Result re = binary(
				function:   ( double E ) => Orbit.meanAnomaly_from_eccentricityAnomaly(E, ecc) - meanAnomaly,
				left:       leftBound,
				right:      rightBound,
				err:        ERROR_THRESHOLD,
				maxIter:    MAX_ITERATION_STEPS
			);
			return re.value;
		}
	}
}
