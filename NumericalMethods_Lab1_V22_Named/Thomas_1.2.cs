public static class Thomas
{
    public static double[] Solve(double[] a, double[] b, double[] c, double[] d)
    {
        int n = b.Length;
        if (n == 0 || a.Length != n || c.Length != n || d.Length != n)
            throw new ArgumentException("Неверные размеры диагоналей");
        if (Math.Abs(a[0]) > 1e-13 || Math.Abs(c[n - 1]) > 1e-13)
            throw new ArgumentException("Граничные диагонали должны быть нулевыми");
        double[] alpha = new double[n], beta = new double[n], x = new double[n];
        for (int i = 0; i < n; i++)
        {
            double denominator = b[i] + (i > 0 ? a[i] * alpha[i - 1] : 0);
            if (Math.Abs(denominator) < 1e-13) throw new InvalidOperationException("Нулевой знаменатель прогонки");
            alpha[i] = i == n - 1 ? 0 : -c[i] / denominator;
            beta[i] = (d[i] - (i > 0 ? a[i] * beta[i - 1] : 0)) / denominator;
        }
        x[n - 1] = beta[n - 1];
        for (int i = n - 2; i >= 0; i--) x[i] = alpha[i] * x[i + 1] + beta[i];
        Console.WriteLine("Прогоночные коэффициенты:");
        for (int i = 0; i < n; i++) Console.WriteLine($"{i + 1}: alpha={alpha[i]:F10} beta={beta[i]:F10}");
        Matrix.Print("x", x);
        double maxResidual = 0;
        for (int i = 0; i < n; i++)
        {
            double sum = b[i] * x[i];
            if (i > 0) sum += a[i] * x[i - 1];
            if (i < n - 1) sum += c[i] * x[i + 1];
            maxResidual = Math.Max(maxResidual, Math.Abs(sum - d[i]));
            Console.WriteLine($"Уравнение {i + 1}: {sum:F10} = {d[i]:F10}");
        }
        Console.WriteLine($"Максимальная невязка = {maxResidual:E6}");
        return x;
    }
}
