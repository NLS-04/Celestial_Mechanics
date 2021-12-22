using System;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Running;
using Vectors;

namespace Celestial_Mechanics {
	public class Rendezvous {
		public static void Main( string[] args ) {
			Orbit a = new Orbit( 1000 * new Vector( -6045, -3490, 2500 ), new Vector( -3457, 6618, 2533 ), 0f, Body.EARTH );
			Orbit b = new Orbit( a.inclination, a.eccentricity, a.semiMajorAxis, a.longitudeOfAscendingNode, a.argumentOfPeriapsis, a.meanAnomaly_At_Epoch, a.epoch, a.body );
			Console.WriteLine( a.staticInformation() );
			Console.WriteLine();
			Console.WriteLine( b.staticInformation() );

			Console.WriteLine(a.Equals(b));

			//Console.WriteLine( Orbit.meanAnomaly( 120, .37255f ) );
			//Console.WriteLine( Orbit.trueAnomaly_to_time( 120, .37255f, 18_834 ) );

			//_ = BenchmarkRunner.Run<Tests>();


			//float maxIntervals = 100f;
			//float ecc = 0f;
			//float me = Constants.pi * 5f/7f;

			//float newton, secant, binary;

			//const int space = -10;

			//string formatStr = "{0," + space + "} | {1," + space + "} | {2," + space + "} | {3," + space + "} | ";

			//Console.WriteLine( String.Format( formatStr, "ecc", "Newton", "Secant", "Binary") );
			//for ( int i = 1 ; i < maxIntervals ; i++ ) {
			//	ecc = i / maxIntervals;

			//	newton = Solver.eccAnomaly_Newton( me, ecc );
			//	secant = Solver.eccAnomaly_Secant( me, ecc );
			//	binary = Solver.eccAnomaly_Binary( me, ecc );

			//	Console.WriteLine( String.Format( formatStr, ecc, newton, secant, binary ) );
			//}
			
			//for ( double t = 0 ; t <= a.period ; t += a.period / 100f ) {
			//	Console.WriteLine( String.Format( "{0, -10} | {1, -10}", t, a.trueAnomaly_at_time(t)*Constants.deg ) );
			//}
		}
	}

	[MemoryDiagnoser]
	public class Tests {
		struct Bassel {
			public int x;
			public int n;
		}
		Bassel bas = new Bassel { x=2, n=5 };


		//[Params(0f, Constants.pi * .5f, Constants.pi * 20f/23f, Constants.pi * 2)]
		[Params( Constants.pi * .5f )]
		public float meanAno;

		[Params( 0f, 0.5f, 0.8f, 1f )]
		//[Params( 0.17f )]
		public float ecc;


		//[Benchmark]
		public void bassel_Jn( ) {
			Solver.bassel_Jn( bas.x, bas.n );
		}

		//[Benchmark]
		public void eccAnomaly_LIT_n5() {
			Solver.LIT_eccAnomaly_n5( meanAno, ecc );
		}

		//[Benchmark]
		public void eccAnomaly_Bassel() {
			Solver.bassel_eccAnomaly( meanAno, ecc );
		}



		//[Benchmark]
		public void eccAnomaly_Bassel_Newton() {
			Solver.eccAnomaly_Bassel_Newton( meanAno, ecc );
		}

		//[Benchmark]
		public void eccAnomaly_LIT_Newton() {
			Solver.eccAnomaly_LIT_Newton( meanAno, ecc );
		}

		[Benchmark]
		public void eccAnomaly_Newton() {
			Solver.eccAnomaly_Newton( meanAno, ecc );
		}

		[Benchmark]
		public void eccAnomaly_Secant() {
			Solver.eccAnomaly_Secant( meanAno, ecc );
		}

		[Benchmark]
		public void eccAnomaly_Binary() {
			Solver.eccAnomaly_Binary( meanAno, ecc );
		}
	}
}