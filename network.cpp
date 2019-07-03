#include<iostream>
#include<vector>
using namespace std;

int main()
{
	int n,m,count;
	cin>>n>>m;
	vector<int> arr(n),a(m),b(m);

	for(int i=0;i<n;i++)
	{cin>>arr[i];}

	for(int j=0;j<m;j++)
	{cin>>a[j]>>b[j];}

	for(int j=0;j<m;j++)
	{
		count=0;
		for(int i=0;i<n;i++)
		{
			if((arr[i]>=a[j])&&(arr[i]<=b[j])){count++;}
		}
		cout<<count<<endl;
	}
	return 0;
}