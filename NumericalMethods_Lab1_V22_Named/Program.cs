using System.Text.Json;

public static class Program
{
    public static void Main()
    {
        Console.OutputEncoding = System.Text.Encoding.UTF8;
        string basePath = AppContext.BaseDirectory;
        using JsonDocument config = JsonDocument.Parse(File.ReadAllText(Path.Combine(basePath, "settings.json")));
        string task = config.RootElement.GetProperty("task").GetString() ?? "all";
        string[] tasks = task == "all" ? ["1.1", "1.2", "1.3", "1.4", "1.5", "1.5_complex"] : [task];
        foreach (string item in tasks)
        {
            Console.WriteLine($"\n================ {item} ================\n");
            try { Run(item, basePath); }
            catch (Exception e) { Console.WriteLine($"ОШИБКА: {e.Message}"); Environment.ExitCode = 1; }
        }
    }

    static void Run(string task, string basePath)
    {
        string name = task.Replace('.', '_');
        string path = Path.Combine(basePath, "Inputs", name + ".json");
        using JsonDocument document = JsonDocument.Parse(File.ReadAllText(path));
        JsonElement root = document.RootElement;
        if (task == "1.2")
        {
            Thomas.Solve(Matrix.Vector(root.GetProperty("a")), Matrix.Vector(root.GetProperty("b")),
                Matrix.Vector(root.GetProperty("c")), Matrix.Vector(root.GetProperty("d")));
            return;
        }
        double[,] a = Matrix.Read(root.GetProperty("A"));
        if (task == "1.1") LU.Run(a, Matrix.Vector(root.GetProperty("b")));
        else if (task == "1.3") Iterations.Run(a, Matrix.Vector(root.GetProperty("b")),
            root.GetProperty("epsilon").GetDouble(), root.GetProperty("maxIterations").GetInt32());
        else if (task == "1.4") JacobiRotation.Run(a,
            root.GetProperty("epsilon").GetDouble(), root.GetProperty("maxIterations").GetInt32());
        else if (task == "1.5" || task == "1.5_complex") QR.Run(a,
            root.GetProperty("epsilon").GetDouble(), root.GetProperty("maxIterations").GetInt32());
        else throw new ArgumentException("Неизвестный номер задания");
    }
}
