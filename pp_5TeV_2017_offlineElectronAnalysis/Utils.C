#include<iostream>
#include<stdio.h>
#include<stdlib.h>
using namespace std;


/**
 * Calculates the greatest common denominator of integers a & b.
 */
int gcd(int a,int b)
{
  if(b==0)
  return a;
  else
  return gcd(b,a%b);
}


/**
 * Prints out the float a as a rational number.
 */
string DoubleRadiansToRationalRadians(float a, bool inUnitsOfPi=true)
{
  if (inUnitsOfPi) a = a/TMath::Pi();
  int t=1;
  while((float)t*a!=(int)(t*a))
  {
   //cout<<t*a<<" "<<(float)t*a<<" "<<(int)(t*a)<<endl; // For Checking
   t=t*10;
  }
  int k=gcd(t*a,t);

  int num = (int)(t*a/k);
  int den = (int)(t/k);
  string numerator;
  string denominator;
  if ((num == 1 && den > 0) || (num == -1 && den < 0)) numerator = "";
  else if ((num == 1 && den < 0) || (num == -1 && den > 0)) numerator = "-";
  else if (num == 0) return "0";
  else numerator = Form("%i", num);

  den = TMath::Abs(den);

  if (den == 1) denominator = "";
  else denominator = Form("/%i", den);
  return numerator + "#pi" + denominator;
}
