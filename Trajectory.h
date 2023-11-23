// Based on u/Rensin2's Desmos solution, translated into UE code by Morwan.
// https://www.reddit.com/r/AskPhysics/comments/17aztlx/trajectory_shifting_calculations_in_space
// Keep in mind this is not a completely finished implementation. 
// The runtime checks will most likely trigger. The results seem to be pretty consistent nevertheless.
// Designed for 2D scenarios in mind.

constexpr int32 NR_MaxIterations = 100;

// https://i.imgur.com/QpdXwks.png
double S_Function(double S, double X, double Y)
{
	const double S_Squared = S * S;
	const double Sqrt_1_Plus_S_Squared = FMath::Sqrt(1 + S_Squared);
	const double C1 = Sqrt_1_Plus_S_Squared * (X + Y * S) + S_Squared;

	return (X * S - Y) * Sqrt_1_Plus_S_Squared - 2 * S - (3 * C1 + 1) / FMath::Sqrt(2 * C1 + 1);
}

// Derivative of the above
double S_FunctionDerivative(double S, double X, double Y)
{
	const double S_Squared = S * S;
	const double Sqrt_1_Plus_S_Squared = FMath::Sqrt(1 + S_Squared);
	const double X_Plus_Y_S = X + Y * S;

	const double C1 = (S * (X * S - Y) / Sqrt_1_Plus_S_Squared);

	const double C2_1 = X_Plus_Y_S * S / Sqrt_1_Plus_S_Squared;
	const double C2_2 = Y * Sqrt_1_Plus_S_Squared + 2 * S;
	const double C2_3 = 3 * (S_Squared + Sqrt_1_Plus_S_Squared * X_Plus_Y_S) + 1;
	const double C2_4 = FMath::Pow(2 * (S_Squared + Sqrt_1_Plus_S_Squared * X_Plus_Y_S) + 1, 1.5);

	const double C2 = ((C2_1 + C2_2) * C2_3) / C2_4;

	const double C3 = X * Sqrt_1_Plus_S_Squared;

	const double C4_1 = 3 * (X_Plus_Y_S * S / Sqrt_1_Plus_S_Squared + Y * Sqrt_1_Plus_S_Squared + 2 * S);
	const double C4_2 = FMath::Sqrt(2 * (S_Squared + Sqrt_1_Plus_S_Squared * X_Plus_Y_S) + 1);

	const double C4 = C4_1 / C4_2;

	return C1 + C2 + C3 - C4 - 2;
}

// Based on https://www.codesansar.com/numerical-methods/newton-raphson-method-using-c-plus-plus.htm
bool S_NewtonRaphsonMethod(double X, double Y, double S_A, double S_B, double& Result)
{
	double S_0 = (S_A + S_B) / 2.0;
	double S_1 = 0.0;
	const double TolerableError = KINDA_SMALL_NUMBER;
	int32 IterationIndex = 0;

	double ClosestSolution = FMath::Abs(S_0);
	double ClosestResult = BIG_NUMBER;

	while (true)
	{
		const double MainFuncResult = S_Function(S_0, X, Y);

		if (MainFuncResult != MainFuncResult || MainFuncResult >= DBL_MAX || MainFuncResult <= -DBL_MAX)
		{
			// Next undefined, return this
			Result = ClosestSolution;
			return true;
		}

		if (FMath::Abs(MainFuncResult) <= TolerableError)
		{
			Result = S_0;
			return true;
		}

		const double MainFuncResultAbs = FMath::Abs(MainFuncResult);
		if (MainFuncResultAbs < ClosestResult)
		{
			ClosestResult = MainFuncResultAbs;
			ClosestSolution = S_0;
		}

		const double DerivativeFuncResult = S_FunctionDerivative(S_0, X, Y);

		if (DerivativeFuncResult == 0.0)
		{
			checkf(false, TEXT("Mathematical error."));
			return false;
		}

		if (MainFuncResult < 0.0)
		{
			S_A = S_0;
		}
		else if (MainFuncResult > 0.0)
		{
			S_B = S_0;
		}

		S_1 = S_0 - MainFuncResult / FMath::Abs(DerivativeFuncResult);

		if (S_1 >= S_A && S_1 <= S_B)
		{
			S_0 = S_1;
		}
		else
		{
			if (MainFuncResult < 0.0)
			{
				S_0 = (S_B + S_0) / 2.0;
			}
			else
			{
				S_0 = (S_A + S_0) / 2.0;
			}
		}

		++IterationIndex;

		if (IterationIndex > NR_MaxIterations)
		{
			Result = ClosestSolution;
			return true;
		}
	}
}

double SMax_Function(double S, double X, double Y)
{
	const double S_Squared = S * S;
	return 2 * (FMath::Sqrt(1 - S_Squared) * (X + Y * S) + S_Squared) + 1;
}

double SMax_FunctionDerivative(double S, double X, double Y)
{
	const double S_Squared = S * S;
	const double Sqrt_One_Minus_S_Squared = FMath::Sqrt(1.0 - S_Squared);

	const double C1 = Y * Sqrt_One_Minus_S_Squared;

	const double C21 = S * (Y * S + X);
	const double C2 = C21 / Sqrt_One_Minus_S_Squared + 2 * S;

	return 2.0 * (C1 - C2);
}

// Based on https://www.codesansar.com/numerical-methods/newton-raphson-method-using-c-plus-plus.htm
bool SMax_NewtonRaphsonMethod(double X, double Y, double& Result)
{
	double S_0 = 0.0;
	double S_1 = 0.0;
	const double TolerableError = KINDA_SMALL_NUMBER;
	int32 IterationIndex = 0;

	double ClosestSolution = FMath::Abs(S_0);
	double ClosestResult = BIG_NUMBER;

	double MainFuncResult = SMax_Function(S_0, X, Y);

	while (true)
	{
		if (FMath::Abs(MainFuncResult) <= TolerableError)
		{
			Result = S_0;
			return true;
		}

		const double DerivativeFuncResult = SMax_FunctionDerivative(S_0, X, Y);

		if (DerivativeFuncResult == 0.0)
		{
			checkf(false, TEXT("Mathematical error."));
			return false;
		}

		S_1 = S_0 - MainFuncResult / DerivativeFuncResult;

		MainFuncResult = SMax_Function(S_1, X, Y);

		if (MainFuncResult != MainFuncResult || MainFuncResult >= DBL_MAX || MainFuncResult <= -DBL_MAX)
		{
			// Next undefined, return this
			Result = ClosestSolution;
			return true;
		}

		const double MainFuncResultAbs = FMath::Abs(MainFuncResult);
		if (MainFuncResultAbs < ClosestResult)
		{
			ClosestResult = MainFuncResultAbs;
			ClosestSolution = S_1;
		}

		S_0 = S_1;

		++IterationIndex;

		if (IterationIndex > NR_MaxIterations)
		{
			Result = ClosestSolution;
			return true;
		}
	}
}

FVector CalculateAutoNavigationAcceleration(const FVector& SpaceshipLocation, const FVector& SpaceshipVelocity, const FVector& TargetLocation, const FVector& TargetVelocity)
{
	FVector Acceleration = FVector::Zero();

	const double A = 9.81;
	const FVector TaregetLocationDiff = TargetLocation - SpaceshipLocation;
	const FVector TargetVelocityDiff = TargetVelocity - SpaceshipVelocity;

	FVector R_N;
	R_N.X = (TaregetLocationDiff.X * TargetVelocityDiff.Y - TaregetLocationDiff.Y * TargetVelocityDiff.X) / TargetVelocityDiff.Length();
	R_N.Y = (TaregetLocationDiff.Y * TargetVelocityDiff.Y + TaregetLocationDiff.X * TargetVelocityDiff.X) / TargetVelocityDiff.Length();

	const double X = (2.0 * A * FMath::Abs(R_N.X)) / FMath::Square(TargetVelocityDiff.Length());
	const double Y = (2.0 * A * R_N.Y) / FMath::Square(TargetVelocityDiff.Length());

	const double S_Min = (Y / X + 1.0 / X);
	double S_Max = 0.0;
	bool bConsiderS_Max = false;

	const double S_A = (Y / X + 1.0 / X);
	const double S_B = (1.0 / X) * (Y + 2 + 2 * FMath::Sqrt(2 + 2 * FMath::Sqrt(FMath::Square(X) + FMath::Square(Y))));

	if (Y < -1)
	{
		if (SMax_NewtonRaphsonMethod(X, Y, S_Max))
		{
			bConsiderS_Max = true;
			S_B = FMath::Min(S_B, S_Max);
		}
		else
		{
			checkf(false, TEXT("Failed to find SMax."));
		}
	}

	double S = 0.0;
	if (S_NewtonRaphsonMethod(X, Y, S_A, S_B, S))
	{
		if (bConsiderS_Max && S > S_Max)
		{
			checkf(false, TEXT("Failed to find S <= SMax."));
		}

		const double One_S_Vector_Length = FVector(1, S, 0).Length();

		const double C_X = TargetVelocityDiff.X * S + (R_N.X <= 0.0 ? -1.0 * TargetVelocityDiff.Y : TargetVelocityDiff.Y);
		const double C_Y = TargetVelocityDiff.Y * S - (R_N.X <= 0.0 ? -1.0 * TargetVelocityDiff.X : TargetVelocityDiff.X);
		const double C_D = One_S_Vector_Length * TargetVelocityDiff.Length();

		Acceleration.X = A * (C_X / C_D);
		Acceleration.Y = A * (C_Y / C_D);
	}

	return Acceleration;
}