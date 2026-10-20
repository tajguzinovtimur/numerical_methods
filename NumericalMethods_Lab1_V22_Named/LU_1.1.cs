public static class LU
{
    public static (double[,] L, double[,] U, double[,] P, int sign) Decompose(double[,] source)
    {
        Matrix.Square(source);
        int n = source.GetLength(0);
        double[,] u = (double[,])source.Clone(), l = Matrix.Identity(n), p = Matrix.Identity(n);
        int sign = 1;
        for (int k = 0; k < n; k++)
        {
            int pivot = k;
            for (int i = k + 1; i < n; i++)
                if (Math.Abs(u[i, k]) > Math.Abs(u[pivot, k])) pivot = i;
            if (Math.Abs(u[pivot, k]) < 1e-13) throw new InvalidOperationException("Нулевой опорный элемент");
            if (pivot != k)
            {
                Swap(u, pivot, k);
                Swap(p, pivot, k);
                for (int j = 0; j < k; j++) (l[pivot, j], l[k, j]) = (l[k, j], l[pivot, j]);
                sign = -sign;
            }
            for (int i = k + 1; i < n; i++)
            {
                l[i, k] = u[i, k] / u[k, k];
                for (int j = k + 1; j < n; j++) u[i, j] -= l[i, k] * u[k, j];
                u[i, k] = 0;
            }
        }
        return (l, u, p, sign);
    }

    static void Swap(double[,] a, int x, int y)
    {
        for (int j = 0; j < a.GetLength(1); j++) (a[x, j], a[y, j]) = (a[y, j], a[x, j]);
    }

    public static double[] Solve(double[,] l, double[,] u, double[,] p, double[] b)
    {
        int n = b.Length;
        double[] pb = Matrix.Multiply(p, b), y = new double[n], x = new double[n];
        for (int i = 0; i < n; i++)
        {
            double sum = pb[i];
            for (int j = 0; j < i; j++) sum -= l[i, j] * y[j];
            y[i] = sum / l[i, i];
        }
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = y[i];
            for (int j = i + 1; j < n; j++) sum -= u[i, j] * x[j];
            x[i] = sum / u[i, i];
        }
        return x;
    }

    public static double[,] Inverse(double[,] l, double[,] u, double[,] p)
    {
        int n = l.GetLength(0);
        double[,] inverse = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            double[] e = new double[n];
            e[j] = 1;
            double[] x = Solve(l, u, p, e);
            for (int i = 0; i < n; i++) inverse[i, j] = x[i];
        }
        return inverse;
    }

    public static void Run(double[,] a, double[] b)
    {
        Matrix.Square(a);
        if (a.GetLength(0) != b.Length) throw new ArgumentException("Неверный размер b");
        var (l, u, p, sign) = Decompose(a);
        double[] x = Solve(l, u, p, b);
        double[,] inverse = Inverse(l, u, p);
        double determinant = sign;
        for (int i = 0; i < b.Length; i++) determinant *= u[i, i];
        Matrix.Print("A", a);
        Matrix.Print("b", b);
        Matrix.Print("L", l);
        Matrix.Print("U", u);
        Matrix.Print("P", p);
        Matrix.Print("L * U", Matrix.Multiply(l, u));
        Matrix.Print("P * A", Matrix.Multiply(p, a));
        Matrix.Print("x", x);
        Matrix.Print("A * x", Matrix.Multiply(a, x));
        Matrix.Print("A^(-1)", inverse);
        Matrix.Print("A * A^(-1)", Matrix.Multiply(a, inverse));
        Console.WriteLine($"\ndet(A) = {determinant:F10}");
        Console.WriteLine($"Ошибка LU - PA = {Matrix.Difference(Matrix.Multiply(l, u), Matrix.Multiply(p, a)):E6}");
        Console.WriteLine($"Ошибка A*A^(-1) - E = {Matrix.Difference(Matrix.Multiply(a, inverse), Matrix.Identity(b.Length)):E6}");
        Console.WriteLine($"Невязка Ax - b = {Matrix.Difference(Matrix.Multiply(a, x), b):E6}");
    }
}
