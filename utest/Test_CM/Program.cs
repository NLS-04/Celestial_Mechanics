using System;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Running;
using Vectors;

namespace Celestial_Mechanics {
	public class Test_CM {
		public static void Main( string[] args ) {
			Orbit a = new Orbit( 1000 * new Vector( -6045, -3490, 2500 ), new Vector( -3457, 6618, 2533 ), 0f, Body.EARTH );
			Orbit b = new Orbit( a.inclination, a.eccentricity, a.semiMajorAxis, a.longitudeOfAscendingNode, a.argumentOfPeriapsis, a.meanAnomaly_At_Epoch, a.epoch, a.body );

			//_ = BenchmarkRunner.Run<Tests_other>();
			//_ = BenchmarkRunner.Run<Tests>();
		}
	}

	[SimpleJob( BenchmarkDotNet.Engines.RunStrategy.Throughput, launchCount: 3, warmupCount: 5, targetCount: 5 )]
	[MemoryDiagnoser]
	public class Tests_other {

		[Params(true, false)]
		public bool searchForMin;

		[Benchmark]
		public void section_search_normal() {
			Solver.section( ( x ) => x / ( x*x + 1 ), -2d, 4d, 36, searchForMin );
		}

		[Benchmark]
		public void section_search_2() {
			Solver.section( ( x ) => x / ( x * x + 1 ), -2d, 4d, 2, searchForMin );
		}

		[Benchmark]
		public void section_search_3() {
			Solver.section( ( x ) => x / ( x * x + 1 ), -2d, 4d, 3, searchForMin );
		}

		[Benchmark]
		public void section_search_4() {
			Solver.section( ( x ) => x / ( x * x + 1 ), -2d, 4d, 4, searchForMin );
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
		public void bassel_Jn() {
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