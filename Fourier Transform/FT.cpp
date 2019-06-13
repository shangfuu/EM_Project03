#include "FT.h"

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // �ť߸��W�v�}�C
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i < M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // �ť߸��W�v�}�C
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag, FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
			//�s�U�ϳť߸���ƻP��Ƴ���
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage,std::complex<double>** FreqRI, int h)
{
	//------------------Initial--------------------------
	int N = h;

	std::vector<std::vector<std::complex<double> >> Freq;
	Freq.resize(N);
	// Freq = InputImage
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			std::complex<double>tmp(InputImage[i][j], 0);
			// std::complex<double>tmp(InputImage[j][i], 0);
			Freq[i].push_back(tmp);
		}
	}

	//-------------------Fast Fourier------------------------

	// Row
	for (int t = 0; t < N; ++t)
	{
		/* bit-reversal permutation */
		for (int i = 1, j = 0; i < N; ++i)
		{
			for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
			if (i > j) swap(Freq[t][i], Freq[t][j]);
		}

		/* dynamic programming */
		for (int k = 2; k <= N; k <<= 1)
		{
			double Omega = -2.0 * PI / k;
			std::complex<double> dSida(cos(Omega), sin(Omega));

			// �Ck�Ӱ��@��FFT
			for (int j = 0; j < N; j += k)
			{
				// �ek/2�ӻP��k/2���T����ƭȫ�n��١A
				// �]������٪��@�_���C
				std::complex<double> Sida(1, 0);
				for (int i = j; i < j + k / 2; i++)
				{
					std::complex<double> a = Freq[t][i];
					std::complex<double> b = Freq[t][i + k / 2] * Sida;
					Freq[t][i] = a + b;
					Freq[t][i + k / 2] = a - b;
					Sida *= dSida;
				}
			}
		}
	}

	// Col
	for (int t = 0; t < N; ++t)
	{
		/* bit-reversal permutation */
		for (int i = 1, j = 0; i < N; ++i)
		{
			for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
			if (i > j) swap(Freq[i][t], Freq[j][t]);
		}

		/* dynamic programming */
		for (int k = 2; k <= N; k <<= 1)
		{
			double Omega = -2.0 * PI / k;
			std::complex<double> dSida(cos(Omega), sin(Omega));

			// �Ck�Ӱ��@��FFT
			for (int j = 0; j < N; j += k)
			{
				// �ek/2�ӻP��k/2���T����ƭȫ�n��١A
				// �]������٪��@�_���C
				std::complex<double> Sida(1, 0);
				for (int i = j; i < j + k / 2; i++)
				{
					std::complex<double> a = Freq[i][t];
					std::complex<double> b = Freq[i + k / 2][t] * Sida;
					Freq[i][t] = a + b;
					Freq[i + k / 2][t] = a - b;
					Sida *= dSida;
				}
			}
		}
	}
	//---------------------Mix Real and Imagine Number----------------------
	int mix;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			FreqRI[i][j].real(Freq[i][j].real() / N);
			FreqRI[i][j].imag(Freq[i][j].imag() / N);

			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			mix = sqrt(pow(FreqRI[i][j].real(), (double) 2.0) + pow(FreqRI[i][j].imag(), (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = mix;
		}
	}

}

std::vector<std::complex<double>> FT::FFT(std::vector<std::complex<double>>x)
{
	// ������ԣ�ofunction��call���N���F
	int N = x.size();

	//------------ Bit-Reversal Permutation--------------------
	for (int i = 1, j = 0; i < N; ++i) {
		for (int k = N >> 1; !((j ^= k)&k); k >>= 1) {
			if (i > j) {
				swap(x[i], x[j]);
			}
		}
	}
	//---------------------Fast Fourier-----------------------
	/* Dynamic Programmin */
	for (int k = 2; k <= N; k <<= 1) {
		double Omega = -2.0 * PI / k;
		std::complex<double> dSida(cos(Omega), sin(Omega));

		// �CK�Ӱ��@��
		for (int j = 0; j < N; j += k) {
			// �e k/2 �ӻP�� k/2 ���T����ƪ���n��١A
			// �]������٤@�_���C
			std::complex<double> Sida(1, 0);
			for (int i = j; i < j + k / 2; i++) {
				std::complex<double>a = x[i];
				std::complex<double>b = x[i + k / 2] * Sida;
				x[i] = a + b;
				x[i + k / 2] = a - b;
				Sida *= dSida;
			}
		}
	}
	return x;
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage,std::complex<double>** FreqRI, int h)
{
	//------------------Initial--------------------------
	int N = h;

	std::vector<std::vector<std::complex<double> >> Freq;
	Freq.resize(N);
	// Freq = preFreq
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Freq[i].push_back(FreqRI[i][j]);
			// Freq[i].push_back(FreqRI[j][i]);
		}
	}

	//-------------------Fast Fourier Inverse------------------------	
	// Row
	for (int t = 0; t < N; ++t)
	{
		/* bit-reversal permutation */
		for (int i = 1, j = 0; i < N; ++i)
		{
			for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
			if (i > j) swap(Freq[t][i], Freq[t][j]);
		}

		/* dynamic programming */
		for (int k = 2; k <= N; k <<= 1)
		{
			double Omega = 2.0 * PI / k;
			std::complex<double> dSida(cos(Omega), sin(Omega));

			// �Ck�Ӱ��@��FFT
			for (int j = 0; j < N; j += k)
			{
				// �ek/2�ӻP��k/2���T����ƭȫ�n��١A
				// �]������٪��@�_���C
				std::complex<double> Sida(1, 0);
				for (int i = j; i < j + k / 2; i++)
				{
					std::complex<double> a = Freq[t][i];
					std::complex<double> b = Freq[t][i + k / 2] * Sida;
					Freq[t][i] = a + b;
					Freq[t][i + k / 2] = a - b;
					Sida *= dSida;
				}
			}
		}
	}

	// Col
	for (int t = 0; t < N; ++t)
	{
		/* bit-reversal permutation */
		for (int i = 1, j = 0; i < N; ++i)
		{
			for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
			if (i > j) swap(Freq[i][t], Freq[j][t]);
		}

		/* dynamic programming */
		for (int k = 2; k <= N; k <<= 1)
		{
			double Omega = 2.0 * PI / k;
			std::complex<double> dSida(cos(Omega), sin(Omega));

			// �Ck�Ӱ��@��FFT
			for (int j = 0; j < N; j += k)
			{
				// �ek/2�ӻP��k/2���T����ƭȫ�n��١A
				// �]������٪��@�_���C
				std::complex<double> Sida(1, 0);
				for (int i = j; i < j + k / 2; i++)
				{
					std::complex<double> a = Freq[i][t];
					std::complex<double> b = Freq[i + k / 2][t] * Sida;
					Freq[i][t] = a + b;
					Freq[i + k / 2][t] = a - b;
					Sida *= dSida;
				}
			}
		}
	}
	//---------------------Mix Real and Imagine Number----------------------
	int mix;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			FreqRI[i][j].real(Freq[i][j].real() / N);
			FreqRI[i][j].imag(Freq[i][j].imag() / N);

			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			mix = sqrt(pow(FreqRI[i][j].real(), (double) 2.0) + pow(FreqRI[i][j].imag(), (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = mix;
		}
	}

}

void FT::InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// �]��FFT �����X�ӡA�ҥH�ڲq�o���Ӥ]�X����
	// ��O�ڴN���������F
}

void FT::LowpassFilter(double** Real, double** Img, int** output, int h, int w)
{
	double **filter = new double*[h];
	for (int i = 0; i < h; i++) {
		filter[i] = new double[w];
	}
	int midh = h / 2;
	int midw = w / 2;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			int u = i - midh, v = j - midw;
			filter[i][j] = 1 / (1 + pow((sqrt(pow(u, 2) + pow(v, 2))) / 10, pow(2, 3)));
			std::complex<double> com(Real[i][j], Img[i][j]);
			com *= filter[i][j];
			Real[i][j] = com.real();
			Img[i][j] = com.imag();
			output[i][j] = sqrt(pow(Real[i][j], 2) + pow(Img[i][j], 2));
		}
	}
}

void FT::HighpassFilter(double** Real, double** Img, int** output, int h, int w)
{
	double **filter = new double*[h];
	for (int i = 0; i < h; i++) {
		filter[i] = new double[w];
	}
	int midh = h / 2;
	int midw = w / 2;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			int u = i - midh, v = j - midw;
			filter[i][j] = 1 - (1 / (1 + pow((sqrt(pow(u, 2) + pow(v, 2))) / 10, pow(2, 3))));
			std::complex<double> com(Real[i][j], Img[i][j]);
			com *= filter[i][j];
			Real[i][j] = com.real();
			Img[i][j] = com.imag();
			output[i][j] = sqrt(pow(Real[i][j], 2) + pow(Img[i][j], 2));
		}
	}
}
